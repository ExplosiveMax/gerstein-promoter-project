import random
random.seed(42)

def gc(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)

def load_fasta(path):
    seqs, cur = [], ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs

# Target GC = mean of both promoter classes combined
allprom = load_fasta("functional_promoters.fasta") + load_fasta("pseudo_promoters.fasta")
region_size = len(allprom[0])
target = sum(gc(s) for s in allprom) / len(allprom)
print(f"Combined promoter mean GC: {target:.3f}, window {region_size}")

lo, hi = target - 0.03, target + 0.03
max_repeat = 0.25
N = 300

print("Loading chr22...")
chr22 = ""
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"):
            chr22 += line.strip()

regions, attempts = [], 0
while len(regions) < N and attempts < 500000:
    attempts += 1
    s = random.randint(0, len(chr22) - region_size)
    seq = chr22[s:s+region_size]
    if sum(1 for c in seq if c.islower())/len(seq) > max_repeat:
        continue
    if not (lo <= gc(seq) <= hi):
        continue
    regions.append((s, s+region_size))

print(f"Found {len(regions)} GC-matched random regions after {attempts} attempts")
with open("random_3way.bed", "w") as out:
    for s, e in regions:
        out.write(f"chr22\t{s}\t{e}\n")
print("Wrote random_3way.bed")