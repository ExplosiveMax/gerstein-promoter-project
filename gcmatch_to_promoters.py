import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tolerance", type=float, default=0.03,
                    help="allowed GC deviation from promoter mean")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

random.seed(args.seed)

def load_fasta(filepath):
    seqs, cur = [], ""
    for line in open(filepath):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper())
            cur = ""
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs

def gc_content(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)

# 1. Read the CURRENT promoters (whatever offset they were extracted at)
promoters = load_fasta("promoters.fasta")
region_size = len(promoters[0])
prom_gc = sum(gc_content(s) for s in promoters) / len(promoters)
print(f"Promoter mean GC: {prom_gc:.3f}, window size: {region_size}")

min_gc = prom_gc - args.tolerance
max_gc = prom_gc + args.tolerance
max_repeat = 0.25
num_regions = 500

# 2. Load chr22
print("Loading chr22...")
chr22_seq = ""
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"):
            chr22_seq += line.strip()

# 3. Sample randoms matched to THIS offset's promoter GC
regions, attempts = [], 0
while len(regions) < num_regions:
    attempts += 1
    start = random.randint(0, len(chr22_seq) - region_size)
    seq = chr22_seq[start:start + region_size]
    if sum(1 for c in seq if c.islower()) / len(seq) > max_repeat:
        continue
    if not (min_gc <= gc_content(seq) <= max_gc):
        continue
    regions.append((start, start + region_size))
    if attempts > 500000:  # safety valve if GC target is too rare
        print("WARNING: hit attempt cap, GC target may be unreachable")
        break

print(f"Found {len(regions)} GC-matched regions after {attempts} attempts")
with open("random_regions.bed", "w") as f:
    for start, end in regions:
        f.write(f"chr22\t{start}\t{end}\n")
print("Saved random_regions.bed")