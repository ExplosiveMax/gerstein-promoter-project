import random

def gc_content(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)

# Load chr22
print("Loading chr22...")
chr22 = ""
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"):
            chr22 += line.strip()

region_size = 1000
num_regions = 500
# target GC ~ your promoter GC, so this mirrors the real experiment
min_gc, max_gc = 0.47, 0.53
max_repeat = 0.25

def sample_regions(rng, n):
    regions = []
    while len(regions) < n:
        start = rng.randint(0, len(chr22) - region_size)
        seq = chr22[start:start + region_size]
        if sum(1 for c in seq if c.islower()) / len(seq) > max_repeat:
            continue
        if not (min_gc <= gc_content(seq) <= max_gc):
            continue
        regions.append((start, start + region_size))
    return regions

# TWO independent random sets, different seeds so they don't overlap
setA = sample_regions(random.Random(111), num_regions)
setB = sample_regions(random.Random(222), num_regions)

# sanity: make sure they don't share coordinates
overlap = set(setA) & set(setB)
print(f"Coordinate overlap between sets: {len(overlap)} (should be 0 or near 0)")

with open("setA.bed", "w") as f:
    for s, e in setA:
        f.write(f"chr22\t{s}\t{e}\n")
with open("setB.bed", "w") as f:
    for s, e in setB:
        f.write(f"chr22\t{s}\t{e}\n")

print("Wrote setA.bed and setB.bed")