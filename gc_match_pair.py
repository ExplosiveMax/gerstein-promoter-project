import random
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cat_a", required=True)
parser.add_argument("--cat_b", required=True)
args = parser.parse_args()
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
a = load_fasta(f"go_{args.cat_a}_pw.fasta")
b = load_fasta(f"go_{args.cat_b}_pw.fasta")

a_bins, b_bins = defaultdict(list), defaultdict(list)
for h, s in a: a_bins[round(gc(s)/BIN)].append((h, s))
for h, s in b: b_bins[round(gc(s)/BIN)].append((h, s))

matched_a, matched_b = [], []
for bin_key in sorted(set(a_bins) | set(b_bins)):
    n = min(len(a_bins[bin_key]), len(b_bins[bin_key]))
    if n == 0: continue
    random.shuffle(a_bins[bin_key]); random.shuffle(b_bins[bin_key])
    matched_a += a_bins[bin_key][:n]
    matched_b += b_bins[bin_key][:n]

def write(path, items):
    with open(path, "w") as f:
        for h, s in items: f.write(f"{h}\n{s}\n")
write(f"go_{args.cat_a}_gcm_{args.cat_b}.fasta", matched_a)
write(f"go_{args.cat_b}_gcm_{args.cat_a}.fasta", matched_b)
gset = sum(gc(s) for _,s in matched_a)/len(matched_a) if matched_a else 0
print(f"{args.cat_a} vs {args.cat_b}: matched {len(matched_a)} each, GC now {gset:.4f}")
