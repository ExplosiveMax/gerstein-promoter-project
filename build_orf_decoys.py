import gzip
import random
from collections import defaultdict
random.seed(42)

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"
UORF_BED = "uorfs_chr22_capped.bed"

STARTS = {"ATG", "CTG", "GTG", "TTG", "ACG"}  # match uORF start diversity
STOPS  = {"TAA", "TAG", "TGA"}

# ---------- load chr22 ----------
print("Loading chr22...")
chr22 = ""
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"):
            chr22 += line.strip()
chr22 = chr22.upper()

def revcomp(s):
    c = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
    return "".join(c.get(b,"N") for b in reversed(s))

# ---------- 5'UTR intervals (strand-aware) ----------
tx = defaultdict(lambda: {"utr": [], "cds": [], "strand": None})
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM or c[2] not in ("UTR","CDS"): continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        tid = None
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("transcript_id"):
                tid = part.split('"')[1]; break
        if tid is None: continue
        tx[tid]["strand"] = strand
        (tx[tid]["utr"] if c[2]=="UTR" else tx[tid]["cds"]).append((start,end))

five = []  # (start,end,strand) 1-based inclusive
for tid,d in tx.items():
    if not d["utr"] or not d["cds"]: continue
    if d["strand"]=="+":
        cds0 = min(s for s,e in d["cds"])
        five += [(s,e,"+") for s,e in d["utr"] if e <= cds0]
    else:
        cds1 = max(e for s,e in d["cds"])
        five += [(s,e,"-") for s,e in d["utr"] if s >= cds1]

# ---------- uORF intervals to exclude ----------
uorf_iv = []
uorf_lengths = []
with open(UORF_BED) as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        s,e = int(c[1]), int(c[2])
        uorf_iv.append((s,e)); uorf_lengths.append(e-s)
uorf_iv.sort()

def overlaps(s,e):
    for us,ue in uorf_iv:
        if us>=e: break
        if ue>s: return True
    return False

# ---------- find decoy ORFs in non-uORF 5'UTR sequence ----------
# For each 5'UTR interval, read its sequence in transcription orientation,
# scan for start->in-frame-stop ORFs, keep those NOT overlapping a real uORF.
decoys = []  # (genomic_start, genomic_end, strand)
for gs, ge, strand in five:
    seq = chr22[gs-1:ge]              # 1-based incl -> slice
    if strand == "-":
        seq = revcomp(seq)
    n = len(seq)
    i = 0
    while i < n-2:
        if seq[i:i+3] in STARTS:
            j = i+3
            while j+3 <= n:
                if seq[j:j+3] in STOPS:
                    orf_len = j+3-i
                    # map local i..j+3 back to genomic coords
                    if strand == "+":
                        g0 = (gs-1) + i
                        g1 = (gs-1) + j+3
                    else:
                        # local coords are on revcomp; map back
                        g1 = ge - i
                        g0 = ge - (j+3)
                    if not overlaps(g0, g1):
                        decoys.append((g0, g1, strand, orf_len))
                    break
                j += 3
        i += 1

print(f"Decoy ORFs found in non-uORF 5'UTR: {len(decoys)}")

# ---------- length-match decoys to the uORF length distribution ----------
# Bucket decoys by length, then for each uORF length pick a decoy of same length.
from collections import defaultdict as dd
by_len = dd(list)
for g0,g1,strand,L in decoys:
    by_len[L].append((g0,g1,strand))

negatives = []
used = set()
fails = 0
for L in uorf_lengths:
    pool = by_len.get(L, [])
    random.shuffle(pool)
    placed = False
    for cand in pool:
        key = (cand[0], cand[1])
        if key in used: continue
        used.add(key); negatives.append(cand); placed=True; break
    if not placed:
        # relax: nearest available length within +-9bp (multiple of 3)
        for dL in [3,-3,6,-6,9,-9]:
            pool = by_len.get(L+dL, [])
            random.shuffle(pool)
            for cand in pool:
                key=(cand[0],cand[1])
                if key in used: continue
                used.add(key); negatives.append(cand); placed=True; break
            if placed: break
    if not placed:
        fails += 1

print(f"Built {len(negatives)} ORF-decoy negatives ({fails} unmatched)")

with open("uorf_decoys_chr22.bed","w") as out:
    for g0,g1,strand in negatives:
        out.write(f"{CHROM}\t{g0}\t{g1}\tdecoy\t0\t{strand}\n")
print("Wrote uorf_decoys_chr22.bed")