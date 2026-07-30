import numpy as np

def load_pfm(path):
    rows = {}
    name = None
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
    # convert counts -> probabilities -> log-odds vs uniform background (0.25 each)
    col_sums = pfm.sum(axis=0)
    probs = (pfm + pseudocount/4) / (col_sums + pseudocount)
    pwm = np.log2(probs / 0.25)
    return pwm

BASE_IDX = {"A":0, "C":1, "G":2, "T":3}
COMP = {"A":"T","T":"A","C":"G","G":"C"}

def revcomp(seq):
    return "".join(COMP[b] for b in reversed(seq))

def score_seq_at(seq_window, pwm):
    L = pwm.shape[1]
    score = 0.0
    for i in range(L):
        b = seq_window[i]
        if b not in BASE_IDX: return -999
        score += pwm[BASE_IDX[b], i]
    return score

def scan(seq, pwm, name):
    L = pwm.shape[1]
    max_score = pwm.max(axis=0).sum()  # theoretical max
    best = (-999, -1, "+")
    for strand, s in [("+", seq), ("-", revcomp(seq))]:
        for i in range(len(s) - L + 1):
            sc = score_seq_at(s[i:i+L], pwm)
            if sc > best[0]:
                pos = i if strand == "+" else len(seq) - L - i
                best = (sc, pos, strand)
    pct = 100 * best[0] / max_score if max_score > 0 else 0
    return best[0], best[1], best[2], pct, max_score

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
    print(f"Loaded {mid} ({name}), length {pwm.shape[1]}")

seqs = load_fasta("onoff_Liver_high.fasta", 5)

print("\n=== Scanning all 5 sequences against all 5 motifs ===")
for si, seq in enumerate(seqs):
    print(f"\n--- Sequence {si} ---")
    for mid, name, pwm in motif_data:
        score, pos, strand, pct, maxscore = scan(seq, pwm, name)
        flag = "  <-- STRONG" if pct > 75 else ("  <-- moderate" if pct > 60 else "")
        print(f"  {name:>8} ({mid}): best score {score:.2f}/{maxscore:.2f} ({pct:.1f}%) at pos {pos} strand {strand}{flag}")
