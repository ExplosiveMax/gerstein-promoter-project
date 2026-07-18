import random
import argparse
from gene_overlap import load_gene_intervals, overlaps_any

parser = argparse.ArgumentParser()
parser.add_argument("--window", type=int, default=500)
parser.add_argument("--gc_tol", type=float, default=0.10)
parser.add_argument("--max_repeat", type=float, default=0.25)
args = parser.parse_args()
random.seed(42)

def gc_content(seq):
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    total = gc + at
    return gc / total if total else 0

def repeat_frac(seq):
    return sum(1 for c in seq if c.islower()) / len(seq)

print("Loading chr22...")
chr22_seq = ""
with open('chr22.fa') as f:
    for line in f:
        if not line.startswith('>'):
            chr22_seq += line.strip()
CHR_LEN = len(chr22_seq)
print(f"chr22 length: {CHR_LEN}")

positives = []
with open('promoters.bed') as f:
    for line in f:
        chrom, s, e, name = line.strip().split('\t')
        positives.append((int(s), int(e)))
print(f"Matching {len(positives)} positives from promoters.bed")

gene_intervals = load_gene_intervals("gencode.v47.basic.annotation.gtf.gz", "chr22")

regions = []
fails = 0
MAX_ATTEMPTS = 1000

for s, e in positives:
    window = e - s
    target_gc = gc_content(chr22_seq[s:e])
    placed = False
    for _ in range(MAX_ATTEMPTS):
        start = random.randint(0, CHR_LEN - window)
        end = start + window
        if overlaps_any(start, end, gene_intervals):
            continue
        seq = chr22_seq[start:end]
        if repeat_frac(seq) > args.max_repeat:
            continue
        if abs(gc_content(seq) - target_gc) <= args.gc_tol:
            regions.append((start, end))
            placed = True
            break
    if not placed:
        fails += 1

print(f"Found {len(regions)} GC-matched (per-sequence, ±{args.gc_tol}), non-genic regions")
print(f"Failed to place: {fails}")
if fails > 0:
    print(f"WARNING: {fails} positives could not get a matched negative — chr22 may be too gene-dense/repeat-heavy at this GC target")

with open('random_regions.bed', 'w') as f:
    for start, end in regions:
        f.write(f'chr22\t{start}\t{end}\n')
print(f"Saved to random_regions.bed with window size {window if positives else args.window}")
