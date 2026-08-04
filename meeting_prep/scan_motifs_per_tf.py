"""Per-TF (individual-motif) ISM-importance overlap breakdown.

Reuses the PFM/PWM/scan/overlap/null functions from
scan_motifs_hivlo_scaleup.py verbatim. The ONLY change is the reporting
granularity: instead of collapsing each panel into one aggregate overlap
rate, we report an overlap rate PER MOTIF (per TF), keeping

  - the liver-specific panel (5 TFs) and the broad panel (14 TFs) SEPARATE
    (never merged into a combined-19 panel), and
  - the high- and low-expression liver sets SEPARATE,

and comparing each motif against its own random-position null.

No model / ISM rerun: this loads the existing ism_importance_{high,low}_seqN.npy
arrays (n=20/class) and the existing onoff_Liver_{high,low}.fasta.

Outputs:
  meeting_prep/motif_per_tf_overlap.csv
  meeting_prep/motif_per_tf_overlap.png
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_SEQS = 20
TOP_N = 15

# --- functions copied unchanged from scan_motifs_hivlo_scaleup.py ----------
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

# --- panels (liver-5 canonical, broad-14 restored) -------------------------
LIVER_MOTIFS = ["MA0060.3", "MA0114.2", "MA0047.4", "MA0102.4", "MA0046.3"]

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

# --- per-motif tally --------------------------------------------------------
def per_motif_rates(panel_name, motif_data):
    """Return list of dicts, one per (motif, expr_set)."""
    # pre-load sequences + importance per expression set
    loaded = {}
    for label, fasta in [("high", "onoff_Liver_high.fasta"),
                         ("low", "onoff_Liver_low.fasta")]:
        seqs = load_fasta(fasta, N_SEQS)
        imps = []
        for si in range(len(seqs)):
            imp_path = f"ism_importance_{label}_seq{si}.npy"
            if not os.path.exists(imp_path):
                raise FileNotFoundError(f"Missing {imp_path}")
            imps.append(np.load(imp_path))
        loaded[label] = (seqs, imps)

    rows = []
    for mid, name, pwm in motif_data:
        for label in ["high", "low"]:
            seqs, imps = loaded[label]
            n_overlap = 0
            null_vals = []
            for si, seq in enumerate(seqs):
                imp = imps[si]
                # deterministic rng, same recipe as scale-up script
                rng = np.random.default_rng(1000 * si + (0 if label == "high" else 1))
                _score, pos, strand, pct, L = scan(seq, pwm)
                overlaps, _ = overlap_with_importance(pos, L, imp, TOP_N)
                if overlaps:
                    n_overlap += 1
                null_vals.append(random_null_rate(L, len(seq), imp, TOP_N, rng=rng))
            rows.append({
                "panel": panel_name,
                "motif_id": mid,
                "tf": name,
                "motif_len": pwm.shape[1],
                "expr_set": label,
                "n_seqs": len(seqs),
                "n_overlap": n_overlap,
                "overlap_rate_pct": round(100 * n_overlap / len(seqs), 1),
                "null_rate_pct": round(100 * float(np.mean(null_vals)), 1),
            })
    return rows


def main():
    liver = load_liver_panel()
    broad = load_broad_panel()
    print(f"liver panel: {len(liver)} TFs   broad panel: {len(broad)} TFs")

    all_rows = per_motif_rates("liver", liver) + per_motif_rates("broad", broad)

    # ---- CSV ----
    import csv
    csv_path = "meeting_prep/motif_per_tf_overlap.csv"
    cols = ["panel", "motif_id", "tf", "motif_len", "expr_set",
            "n_seqs", "n_overlap", "overlap_rate_pct", "null_rate_pct"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"wrote {csv_path}  ({len(all_rows)} rows)")

    # ---- console table ----
    print(f"\n{'panel':>6} {'TF':>10} {'motif':>9} {'set':>5} "
          f"{'overlap%':>9} {'null%':>7}")
    for r in all_rows:
        print(f"{r['panel']:>6} {r['tf']:>10} {r['motif_id']:>9} "
              f"{r['expr_set']:>5} {r['overlap_rate_pct']:>8.1f}% "
              f"{r['null_rate_pct']:>6.1f}%")

    # ---- bar chart: separate subplot per panel; per-TF high/low bars + null marker ----
    C_HIGH, C_LOW, C_NULL = "#2b6cb0", "#dd6b20", "#4a5568"
    for panel, motif_data in [("liver", liver), ("broad", broad)]:
        pass  # panels drawn below together

    fig, axes = plt.subplots(1, 2, figsize=(16, 6),
                             gridspec_kw={"width_ratios": [len(liver), len(broad)]})
    for ax, panel in zip(axes, ["liver", "broad"]):
        prows = [r for r in all_rows if r["panel"] == panel]
        tfs = sorted({r["tf"] for r in prows})
        # keep original panel order
        order = LIVER_MOTIFS if panel == "liver" else \
            sorted({r["motif_id"] for r in prows})
        id2tf = {r["motif_id"]: r["tf"] for r in prows}
        tfs = [id2tf[m] for m in order]
        x = np.arange(len(tfs))
        w = 0.38
        def get(mid, label, key):
            return next(r[key] for r in prows
                       if r["motif_id"] == mid and r["expr_set"] == label)
        high = [get(m, "high", "overlap_rate_pct") for m in order]
        low = [get(m, "low", "overlap_rate_pct") for m in order]
        null = [(get(m, "high", "null_rate_pct") + get(m, "low", "null_rate_pct")) / 2
                for m in order]
        ax.bar(x - w/2, high, w, label="high overlap%", color=C_HIGH)
        ax.bar(x + w/2, low, w, label="low overlap%", color=C_LOW)
        # null as a horizontal tick spanning each TF group
        for xi, nv in zip(x, null):
            ax.plot([xi - w, xi + w], [nv, nv], color=C_NULL, lw=2,
                    solid_capstyle="butt",
                    label="random-position null" if xi == 0 else None)
        ax.set_xticks(x)
        ax.set_xticklabels(tfs, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("overlap with top-15 ISM positions (%)")
        ax.set_title(f"{panel}-specific panel ({len(tfs)} TFs)" if panel == "liver"
                     else f"broad TF panel ({len(tfs)} TFs)")
        ax.set_ylim(0, 100)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("Per-TF ISM-importance overlap vs random-position null "
                 "(liver high vs low, n=20/class, top-15 positions)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png = "meeting_prep/motif_per_tf_overlap.png"
    fig.savefig(png, dpi=300)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
