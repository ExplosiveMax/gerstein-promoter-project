import gzip
import random

random.seed(42)

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"
WINDOW = 500
N = 500

# 1. Collect gene bodies on chr22 (start, end), filtering to reasonably long ones
genes = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene" or c[0] != CHROM:
            continue
        start, end = int(c[3]), int(c[4])
        if end - start >= WINDOW + 2000:  # long enough to sample an interior window
            genes.append((start, end))

print(f"Genes on {CHROM} long enough to sample: {len(genes)}")

# 2. Load chr22 to check repeat content of sampled windows
print("Loading chr22...")
chr22 = ""
with open("chr22.fa") as fh:
    for line in fh:
        if not line.startswith(">"):
            chr22 += line.strip()

# 3. Sample interior windows: pick a random gene, random start well inside it.
#    Avoid the first/last 1000bp so we're NOT sitting on the promoter/TSS or the very end.
regions = []
attempts = 0
while len(regions) < N and attempts < 200000:
    attempts += 1
    gstart, gend = random.choice(genes)
    lo = gstart + 1000
    hi = gend - 1000 - WINDOW
    if hi <= lo:
        continue
    s = random.randint(lo, hi)
    seq = chr22[s:s+WINDOW]
    # skip mostly-repetitive windows, same threshold as elsewhere
    if sum(1 for ch in seq if ch.islower()) / len(seq) > 0.25:
        continue
    regions.append((s, s+WINDOW))

print(f"Sampled {len(regions)} interior genic windows after {attempts} attempts")

with open("genic_windows.bed", "w") as out:
    for s, e in regions:
        out.write(f"{CHROM}\t{s}\t{e}\n")
print("Wrote genic_windows.bed")