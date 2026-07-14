import gzip
from collections import defaultdict

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"

# --- arbitrary choices you may want to revisit with Joel ---
MIN_UORF_LEN = 6      # min length in codons (incl. start+stop). 6 codons = 18nt. Some papers use different floors.
START_CODONS = {"ATG"}  # canonical only. Some uORF studies include near-cognate starts (CTG, GTG, etc.)
STOP_CODONS = {"TAA", "TAG", "TGA"}
# -----------------------------------------------------------

# 1. Collect UTR, CDS, exon features per transcript on chr22
tx = defaultdict(lambda: {"utr": [], "cds": [], "strand": None})

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM:
            continue
        feat = c[2]
        if feat not in ("UTR", "CDS"):
            continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        # pull transcript_id
        tid = None
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("transcript_id"):
                tid = part.split('"')[1]
                break
        if tid is None:
            continue
        tx[tid]["strand"] = strand
        if feat == "UTR":
            tx[tid]["utr"].append((start, end))
        else:
            tx[tid]["cds"].append((start, end))

print(f"Transcripts on {CHROM} with UTR+CDS: "
      f"{sum(1 for t in tx.values() if t['utr'] and t['cds'])}")

# 2. For each transcript, figure out which UTR pieces are 5' (before CDS start)
#    On + strand: 5'UTR is UTR with end <= min(CDS start)
#    On - strand: 5'UTR is UTR with start >= max(CDS end)
five_utrs = {}  # tid -> list of (start,end) intervals that are 5'UTR

for tid, d in tx.items():
    if not d["utr"] or not d["cds"]:
        continue
    strand = d["strand"]
    if strand == "+":
        cds_start = min(s for s, e in d["cds"])
        utr5 = [(s, e) for s, e in d["utr"] if e <= cds_start]
    else:
        cds_end = max(e for s, e in d["cds"])
        utr5 = [(s, e) for s, e in d["utr"] if s >= cds_end]
    if utr5:
        five_utrs[tid] = (strand, sorted(utr5))

print(f"Transcripts with an identifiable 5'UTR: {len(five_utrs)}")

# 3. Load chr22 sequence
print("Loading chr22...")
seq = []
with open("chr22.fa") as f:
    for line in f:
        if not line.startswith(">"):
            seq.append(line.strip())
chr22 = "".join(seq)

def revcomp(s):
    comp = {"A":"T","T":"A","C":"G","G":"C","N":"N",
            "a":"t","t":"a","c":"g","g":"c","n":"n"}
    return "".join(comp.get(b, "N") for b in reversed(s))

# 4. Build the 5'UTR sequence (in transcription direction) and scan for ORFs
def get_utr5_seq(strand, intervals):
    # intervals are genomic; BED-like half-open not assumed here (GTF is 1-based inclusive)
    pieces = [chr22[s-1:e] for s, e in intervals]  # 1-based inclusive -> python slice
    s = "".join(pieces)
    if strand == "-":
        s = revcomp(s)
    return s.upper()

def find_uorfs(utr_seq):
    orfs = []
    n = len(utr_seq)
    for i in range(n - 2):
        if utr_seq[i:i+3] in START_CODONS:
            # scan in frame for a stop
            j = i + 3
            while j + 3 <= n:
                codon = utr_seq[j:j+3]
                if codon in STOP_CODONS:
                    length_codons = (j + 3 - i) // 3
                    if length_codons >= MIN_UORF_LEN:
                        orfs.append(utr_seq[i:j+3])
                    break
                j += 3
    return orfs

all_uorfs = []
for tid, (strand, intervals) in five_utrs.items():
    utr_seq = get_utr5_seq(strand, intervals)
    if "N" in utr_seq or len(utr_seq) < MIN_UORF_LEN * 3:
        continue
    for orf in find_uorfs(utr_seq):
        all_uorfs.append((tid, orf))

print(f"Candidate uORFs found: {len(all_uorfs)}")

# 5. Write to FASTA (dedup identical sequences)
seen = set()
written = 0
with open("uorfs.fasta", "w") as out:
    for tid, orf in all_uorfs:
        if orf in seen:
            continue
        seen.add(orf)
        out.write(f">{tid}_uorf_{written}\n{orf}\n")
        written += 1

print(f"Wrote {written} unique uORF sequences to uorfs.fasta")