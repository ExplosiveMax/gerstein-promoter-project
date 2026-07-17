import gzip, random
from collections import defaultdict
random.seed(42)

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"

def load_lens(f):
    seqs, cur = [], ""
    for l in open(f):
        l = l.strip()
        if l.startswith(">"):
            if cur: seqs.append(len(cur)); cur = ""
        else: cur += l
    if cur: seqs.append(len(cur))
    return seqs

target_lens = load_lens("uorfs_capped.fasta")
print(f"Matching {len(target_lens)} uORF lengths")

print("Loading chr22...")
chr22 = ""
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"): chr22 += line.strip()

# collect UTR and CDS intervals 
utr_iv, cds_iv = [], []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM: continue
        if c[2] == "UTR": utr_iv.append((int(c[3])-1, int(c[4])))
        elif c[2] == "CDS": cds_iv.append((int(c[3])-1, int(c[4])))

utr_iv.sort()
print(f"UTR intervals: {len(utr_iv)}, CDS intervals: {len(cds_iv)}")

def hits_utr(s, e):
    for us, ue in utr_iv:
        if us >= e: break
        if ue > s: return True
    return False

def repeat_frac(seq):
    return sum(1 for c in seq if c.islower()) / len(seq)

# CLASS 3: non-UTR (CDS sequence, no UTR overlap) 
cds_pool = [(s, e) for s, e in cds_iv if e - s >= 300]
nonutr = []
fails = 0
for L in target_lens:
    placed = False
    for _ in range(200):
        s0, e0 = random.choice(cds_pool)
        if e0 - L <= s0: continue
        ws = random.randint(s0, e0 - L); we = ws + L
        if hits_utr(ws, we): continue
        if repeat_frac(chr22[ws:we]) > 0.25: continue
        nonutr.append((ws, we)); placed = True; break
    if not placed: fails += 1
print(f"non-UTR: {len(nonutr)} built ({fails} failed)")

# --- CLASS 4: random (length-matched, repeat-filtered) ---
rand = []
fails = 0
for L in target_lens:
    placed = False
    for _ in range(200):
        ws = random.randint(0, len(chr22) - L); we = ws + L
        if repeat_frac(chr22[ws:we]) > 0.25: continue
        rand.append((ws, we)); placed = True; break
    if not placed: fails += 1
print(f"random: {len(rand)} built ({fails} failed)")

for name, regs in [("nonutr_4way", nonutr), ("random_4way", rand)]:
    with open(f"{name}.bed", "w") as out:
        for s, e in regs: out.write(f"{CHROM}\t{s}\t{e}\t{name}\t0\t+\n")
    print(f"Wrote {name}.bed")