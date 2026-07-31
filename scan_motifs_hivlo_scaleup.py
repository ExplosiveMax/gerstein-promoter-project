"""Scaled-up motif / ISM-importance overlap scan.

Adapted from scan_motifs_hivlo.py (n=5, 5 liver motifs). Changes:
  - n=20 high + 20 low sequences (set via --n_seqs)
  - Two motif panels: the original 5 liver-specific TFs, and a broader
    panel of ~14 broadly-acting human TFs (jaspar_motifs_broad/), plus a
    combined run, to rule out "the model uses different motifs than the 5
    liver ones I picked."
  - Adds a random-position null: for each (seq, motif) pair we also drop the
    motif window at a random location and check overlap, averaged over many
    trials, to report the chance overlap rate for the same window sizes.

This is a scale-up of a NEGATIVE result, not a search for a positive one:
we run the SAME trained checkpoint (onoff_liver_model_seed42) and the SAME
importance arrays, only widening n and the motif panel.
"""
import argparse
import glob
import os
import numpy as np

# ---------------------------------------------------------------------------
# PFM / PWM helpers (unchanged from scan_motifs_hivlo.py)
# ---------------------------------------------------------------------------
def load_pfm(path):
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line.split("\t")[1] if "\t" in line else line[1:]
            elif line:
                base = line[0]
                nums = line[line.index("[") + 1 : line.index("]")].split()
                rows[base] = [int(x) for x in nums]
    pfm = np.array([rows["A"], rows["C"], rows["G"], rows["T"]], dtype=float)
    return name, pfm

def pfm_to_pwm(pfm, pseudocount=0.8):
    col_sums = pfm.sum(axis=0)
    probs = (pfm + pseudocount / 4) / (col_sums + pseudocount)
    return np.log2(probs / 0.25)

BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def revcomp(seq):
    return "".join(COMP[b] for b in reversed(seq))

def score_at(window, pwm):
    L = pwm.shape[1]
    s = 0.0
    for i in range(L):
        b = window[i]
        if b not in BASE_IDX:
            return -999
        s += pwm[BASE_IDX[b], i]
    return s

def scan(seq, pwm):
    L = pwm.shape[1]
    max_score = pwm.max(axis=0).sum()
    best = (-999, -1, "+")
    for strand, s in [("+", seq), ("-", revcomp(seq))]:
        for i in range(len(s) - L + 1):
            sc = score_at(s[i : i + L], pwm)
            if sc > best[0]:
                pos = i if strand == "+" else len(seq) - L - i
                best = (sc, pos, strand)
    pct = 100 * best[0] / max_score if max_score > 0 else 0
    return best[0], best[1], best[2], pct, L

def load_fasta(filepath, n):
    seqs, cur = [], ""
    for line in open(filepath):
        line = line.strip()
        if line.startswith(">"):
            if cur:
                seqs.append(cur.upper())
                cur = ""
        else:
            cur += line
    if cur:
        seqs.append(cur.upper())
    return seqs[:n]

def overlap_with_importance(pos, motif_len, importance, top_n=15):
    top_positions = set(np.argsort(importance)[::-1][:top_n])
    motif_positions = set(range(max(0, pos), pos + motif_len))
    inter = top_positions & motif_positions
    return len(inter) > 0, len(inter)

def random_null_rate(motif_len, seq_len, importance, top_n=15, trials=2000, rng=None):
    """Chance probability that a motif window of this length, placed at a
    uniformly random position, overlaps the top_n importance positions."""
    if rng is None:
        rng = np.random.default_rng(0)
    top_positions = set(np.argsort(importance)[::-1][:top_n])
    hits = 0
    max_start = seq_len - motif_len
    for _ in range(trials):
        pos = int(rng.integers(0, max_start + 1))
        if top_positions & set(range(pos, pos + motif_len)):
            hits += 1
    return hits / trials

# ---------------------------------------------------------------------------
# Motif panels
# ---------------------------------------------------------------------------
LIVER_MOTIFS = ["MA0060.3", "MA0114.2", "MA0047.4", "MA0102.4", "MA0046.3"]  # NFYA, HNF4A, FOXA2, CEBPA, HNF1A

def load_liver_panel():
    out = []
    for mid in LIVER_MOTIFS:
        name, pfm = load_pfm(f"jaspar_motifs/{mid}.pfm")
        out.append((mid, name, pfm_to_pwm(pfm)))
    return out

def load_broad_panel():
    out = []
    for path in sorted(glob.glob("jaspar_motifs_broad/*.pfm")):
        mid = os.path.splitext(os.path.basename(path))[0]
        name, pfm = load_pfm(path)
        out.append((mid, name, pfm_to_pwm(pfm)))
    return out

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run_panel(panel_name, motif_data, n_seqs, top_n, verbose=False):
    results = []
    for label, fasta in [("high", "onoff_Liver_high.fasta"), ("low", "onoff_Liver_low.fasta")]:
        seqs = load_fasta(fasta, n_seqs)
        for si, seq in enumerate(seqs):
            imp_path = f"ism_importance_{label}_seq{si}.npy"
            if not os.path.exists(imp_path):
                raise FileNotFoundError(
                    f"Missing {imp_path}. Run run_ism.py with --n_seqs {n_seqs} first."
                )
            imp = np.load(imp_path)
            rng = np.random.default_rng(1000 * si + (0 if label == "high" else 1))
            for mid, name, pwm in motif_data:
                _score, pos, strand, pct, L = scan(seq, pwm)
                overlaps, n_overlap = overlap_with_importance(pos, L, imp, top_n)
                null = random_null_rate(L, len(seq), imp, top_n, rng=rng)
                results.append((label, si, name, mid, pct, pos, strand, overlaps, n_overlap, null))

    print(f"\n{'='*70}\nPANEL: {panel_name}   ({len(motif_data)} motifs, "
          f"top_n={top_n} importance positions, n_seqs={n_seqs}/class)\n{'='*70}")
    if verbose:
        print(f"{'label':>5} {'seq':>3} {'motif':>8} {'match%':>7} {'pos':>5} {'strand':>6} {'overlaps':>18}")
        for r in sorted(results, key=lambda x: (-x[4])):
            label, si, name, mid, pct, pos, strand, overlaps, n_overlap, null = r
            flag = f"YES ({n_overlap} pos)" if overlaps else "no"
            print(f"{label:>5} {si:>3} {name:>8} {pct:>6.1f}% {pos:>5} {strand:>6} {flag:>18}")

    summary = {}
    for label in ["high", "low"]:
        subset = [r for r in results if r[0] == label]
        n_ov = sum(1 for r in subset if r[7])
        mean_null = float(np.mean([r[9] for r in subset]))
        rate = 100 * n_ov / len(subset)
        summary[label] = (n_ov, len(subset), rate, 100 * mean_null)
        print(f"  {label}: {n_ov}/{len(subset)} motif-hits overlap top-{top_n} "
              f"importance ({rate:.1f}%)   |  random-position null: {100*mean_null:.1f}%")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_seqs", type=int, default=20)
    ap.add_argument("--top_n", type=int, default=15)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    liver = load_liver_panel()
    broad = load_broad_panel()
    combined = liver + broad

    s_liver = run_panel("liver-specific (original 5)", liver, args.n_seqs, args.top_n)
    s_broad = run_panel("broad TF panel (added)", broad, args.n_seqs, args.top_n)
    s_comb = run_panel("combined (liver + broad)", combined, args.n_seqs, args.top_n)

    print(f"\n{'#'*70}\nOVERALL SUMMARY (n_seqs={args.n_seqs}/class, top_n={args.top_n})\n{'#'*70}")
    for pname, s in [("liver-5", s_liver), ("broad-14", s_broad), ("combined-19", s_comb)]:
        for label in ["high", "low"]:
            n_ov, tot, rate, null = s[label]
            print(f"  {pname:>12} | {label:>4}: {n_ov:>3}/{tot:<3} = {rate:5.1f}%   (null {null:.1f}%)")


if __name__ == "__main__":
    main()
