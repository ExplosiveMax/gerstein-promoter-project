import numpy as np

def load_pfm(path):
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line.split("\t")[1] if "\t" in line else line[1:]
            elif line:
                base = line[0]
                nums = line[line.index("[")+1 : line.index("]")].split()
                rows[base] = [int(x) for x in nums]
    pfm = np.array([rows["A"], rows["C"], rows["G"], rows["T"]], dtype=float)
    return name, pfm

def pfm_to_pwm(pfm, pseudocount=0.8):
    col_sums = pfm.sum(axis=0)
    probs = (pfm + pseudocount/4) / (col_sums + pseudocount)
    return np.log2(probs / 0.25)

BASE_IDX = {"A":0,"C":1,"G":2,"T":3}
COMP = {"A":"T","T":"A","C":"G","G":"C"}
def revcomp(seq): return "".join(COMP[b] for b in reversed(seq))

def score_at(window, pwm):
    L = pwm.shape[1]; s = 0.0
    for i in range(L):
        b = window[i]
        if b not in BASE_IDX: return -999
        s += pwm[BASE_IDX[b], i]
    return s

def scan(seq, pwm):
    L = pwm.shape[1]
    max_score = pwm.max(axis=0).sum()
    best = (-999, -1, "+")
    for strand, s in [("+", seq), ("-", revcomp(seq))]:
        for i in range(len(s) - L + 1):
            sc = score_at(s[i:i+L], pwm)
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
            if cur: seqs.append(cur.upper()); cur = ""
        else: cur += line
    if cur: seqs.append(cur.upper())
    return seqs[:n]

MOTIFS = ["MA0060.3", "MA0114.2", "MA0047.4", "MA0102.4", "MA0046.3"]
motif_data = []
for mid in MOTIFS:
    name, pfm = load_pfm(f"jaspar_motifs/{mid}.pfm")
    pwm = pfm_to_pwm(pfm)
    motif_data.append((mid, name, pwm))

def overlap_with_importance(pos, motif_len, importance, top_n=15):
    top_positions = set(np.argsort(importance)[::-1][:top_n])
    motif_positions = set(range(max(0,pos), pos+motif_len))
    return len(top_positions & motif_positions) > 0, len(top_positions & motif_positions)

results = []
for label, fasta in [("high", "onoff_Liver_high.fasta"), ("low", "onoff_Liver_low.fasta")]:
    seqs = load_fasta(fasta, 5)
    for si, seq in enumerate(seqs):
        imp = np.load(f"ism_importance_{label}_seq{si}.npy")
        for mid, name, pwm in motif_data:
            score, pos, strand, pct, L = scan(seq, pwm)
            overlaps, n_overlap = overlap_with_importance(pos, L, imp)
            results.append((label, si, name, mid, pct, pos, strand, overlaps, n_overlap))

print(f"{'label':>5} {'seq':>3} {'motif':>8} {'match%':>7} {'pos':>5} {'strand':>6} {'overlaps_top15_importance':>26}")
for r in sorted(results, key=lambda x: (-x[4])):
    label, si, name, mid, pct, pos, strand, overlaps, n_overlap = r
    flag = f"YES ({n_overlap} pos)" if overlaps else "no"
    print(f"{label:>5} {si:>3} {name:>8} {pct:>6.1f}% {pos:>5} {strand:>6} {flag:>26}")

print("\n=== Summary: overlap rate by expression level ===")
for label in ["high", "low"]:
    subset = [r for r in results if r[0] == label]
    n_overlap = sum(1 for r in subset if r[7])
    print(f"{label}: {n_overlap}/{len(subset)} motif-hits overlap top importance positions ({100*n_overlap/len(subset):.1f}%)")
