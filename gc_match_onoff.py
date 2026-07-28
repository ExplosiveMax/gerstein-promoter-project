import random
from collections import defaultdict
random.seed(42)

def load_fasta(path):
    out, cur, hdr = [], "", None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: out.append((hdr, cur.upper()))
            hdr, cur = line, ""
        else: cur += line
    if cur: out.append((hdr, cur.upper()))
    return out

def gc(seq):
    g = seq.count("G")+seq.count("C"); at = seq.count("A")+seq.count("T")
    return g/(g+at) if (g+at) else 0

BIN = 0.02
for tissue in ["Liver", "Muscle_Skeletal", "Whole_Blood"]:
    high = load_fasta(f"onoff_{tissue}_high.fasta")
    low  = load_fasta(f"onoff_{tissue}_low.fasta")

    high_bins, low_bins = defaultdict(list), defaultdict(list)
    for h, s in high: high_bins[round(gc(s)/BIN)].append((h, s))
    for h, s in low:  low_bins[round(gc(s)/BIN)].append((h, s))

    matched_high, matched_low = [], []
    for b in sorted(set(high_bins) | set(low_bins)):
        n = min(len(high_bins[b]), len(low_bins[b]))
        if n == 0: continue
        random.shuffle(high_bins[b]); random.shuffle(low_bins[b])
        matched_high += high_bins[b][:n]
        matched_low += low_bins[b][:n]

    def write(path, items):
        with open(path, "w") as f:
            for h, s in items: f.write(f"{h}\n{s}\n")
    write(f"onoff_{tissue}_high_gcmatched.fasta", matched_high)
    write(f"onoff_{tissue}_low_gcmatched.fasta", matched_low)

    hg = sum(gc(s) for _,s in matched_high)/len(matched_high) if matched_high else 0
    lg = sum(gc(s) for _,s in matched_low)/len(matched_low) if matched_low else 0
    print(f"{tissue}: matched {len(matched_high)} high, {len(matched_low)} low")
    print(f"  mean GC after matching: high={hg:.4f} low={lg:.4f}")
