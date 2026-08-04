import gzip, random, argparse
from pyfaidx import Fasta

p = argparse.ArgumentParser()
p.add_argument("--offset", type=int, default=0)
p.add_argument("--window", type=int, default=500)
p.add_argument("--limit", type=int, default=500)
p.add_argument("--seed", type=int, default=42)
p.add_argument("--genome", default="/Users/maxwellkim/promoter_project/.claude/worktrees/awesome-mendeleev-f66987/GRCh38.primary_assembly.genome.fa")
args = p.parse_args()
random.seed(args.seed)

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr1"
HALF = args.window // 2
MAXOFF = 1000000          # reserve room so the SAME genes work at every offset
GC_TOL = 0.10
MAX_REPEAT = 0.25
MAX_N = 0.05

genome = Fasta(args.genome)
CHR_LEN = len(genome[CHROM])

def gc(s):
    s = s.upper(); g = s.count("G")+s.count("C"); a = s.count("A")+s.count("T")
    return g/(g+a) if (g+a) else 0
def frac_lower(s): return sum(1 for c in s if c.islower())/len(s)
def frac_n(s): return s.upper().count("N")/len(s)

# ---- fixed gene set: usable at EVERY offset, no isolation constraint ----
genes = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene" or c[0] != CHROM: continue
        s, e, strand = int(c[3]), int(c[4]), c[6]
        tss = s if strand == "+" else e
        lo = tss - MAXOFF - HALF
        hi = tss + MAXOFF + HALF
        if lo < 0 or hi > CHR_LEN: continue
        genes.append((tss, strand))
random.shuffle(genes)
genes = genes[:args.limit]
print(f"Fixed gene set: {len(genes)} genes on {CHROM}")

# ---- positives at this offset ----
pos = []
for tss, strand in genes:
    center = tss - args.offset if strand == "+" else tss + args.offset
    ws = center - HALF; we = ws + args.window
    seq = str(genome[CHROM][ws:we])
    if frac_n(seq) > MAX_N: continue
    pos.append(seq)
print(f"Positives: {len(pos)}")

# ---- GC-matched negatives, sampled anywhere on chr1 ----
neg, fails = [], 0
for s in pos:
    target = gc(s); placed = False
    for _ in range(3000):
        ws = random.randint(0, CHR_LEN - args.window); we = ws + args.window
        cand = str(genome[CHROM][ws:we])
        if frac_n(cand) > MAX_N: continue
        if frac_lower(cand) > MAX_REPEAT: continue
        if abs(gc(cand) - target) <= GC_TOL:
            neg.append(cand); placed = True; break
    if not placed: fails += 1
print(f"Negatives: {len(neg)} ({fails} unmatched)")

with open("promoters.fasta","w") as f:
    for i,s in enumerate(pos): f.write(f">pos_{i}\n{s.upper()}\n")
with open("random_sequences.fa","w") as f:
    for i,s in enumerate(neg): f.write(f">neg_{i}\n{s.upper()}\n")
print("Wrote promoters.fasta / random_sequences.fa")
