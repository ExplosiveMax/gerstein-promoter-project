import gzip
import random
random.seed(42)

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"
UORF_BED = "uorfs_chr22_capped.bed"

# ---------- 1. Collect 5'UTR intervals on chr22 from the GTF ----------
# 5'UTR = UTR features upstream of the CDS start (strand-aware)
from collections import defaultdict
tx = defaultdict(lambda: {"utr": [], "cds": [], "strand": None})

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM or c[2] not in ("UTR", "CDS"):
            continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        tid = None
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("transcript_id"):
                tid = part.split('"')[1]; break
        if tid is None:
            continue
        tx[tid]["strand"] = strand
        (tx[tid]["utr"] if c[2] == "UTR" else tx[tid]["cds"]).append((start, end))

five_utr = []  # (start, end) genomic intervals of 5'UTRs, 1-based inclusive
for tid, d in tx.items():
    if not d["utr"] or not d["cds"]:
        continue
    if d["strand"] == "+":
        cds_start = min(s for s, e in d["cds"])
        five_utr += [(s, e) for s, e in d["utr"] if e <= cds_start]
    else:
        cds_end = max(e for s, e in d["cds"])
        five_utr += [(s, e) for s, e in d["utr"] if s >= cds_end]

print(f"5'UTR intervals on {CHROM}: {len(five_utr)}")

# ---------- 2. Load uORF intervals (to exclude) and their lengths ----------
uorf_intervals = []
uorf_lengths = []
with open(UORF_BED) as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        s, e = int(c[1]), int(c[2])          # BED is 0-based half-open
        uorf_intervals.append((s, e))
        uorf_lengths.append(e - s)
print(f"uORFs to match: {len(uorf_lengths)}")

# ---------- 3. Build a set of genomic positions covered by uORFs ----------
# (for overlap checking). Use interval list; check overlap by scanning.
uorf_intervals.sort()

def overlaps_uorf(s, e):
    # simple overlap check against sorted uORF intervals
    for us, ue in uorf_intervals:
        if us >= e:      # past the window, no more overlaps possible
            break
        if ue > s:       # us < e already implied; ue>s means overlap
            return True
    return False

# ---------- 4. Build pool of non-uORF 5'UTR sub-intervals ----------
# Convert 5'UTR intervals (1-based incl) to 0-based half-open, then
# keep only the parts not overlapping any uORF.
def subtract_uorfs(s, e):
    """Return list of (start,end) sub-intervals of [s,e) not covered by uORFs."""
    pieces = [(s, e)]
    for us, ue in uorf_intervals:
        if ue <= s or us >= e:
            continue
        new = []
        for ps, pe in pieces:
            if ue <= ps or us >= pe:
                new.append((ps, pe))
            else:
                if ps < us: new.append((ps, us))
                if ue < pe: new.append((ue, pe))
        pieces = new
    return pieces

non_uorf = []
for s, e in five_utr:
    s0 = s - 1  # to 0-based half-open
    for ps, pe in subtract_uorfs(s0, e):
        if pe - ps >= 12:   # keep pieces at least min uORF length
            non_uorf.append((ps, pe))

total_bp = sum(pe - ps for ps, pe in non_uorf)
print(f"Non-uORF 5'UTR pieces: {len(non_uorf)}, total {total_bp} bp")

# ---------- 5. For each uORF length, sample a length-matched window ----------
# Sort pieces by length so we can find ones big enough.
non_uorf.sort(key=lambda x: x[1]-x[0], reverse=True)
random.shuffle(non_uorf)

negatives = []
used = set()  # avoid reusing exact same window
fails = 0
for L in uorf_lengths:
    placed = False
    # try random pieces that are long enough
    candidates = [p for p in non_uorf if (p[1]-p[0]) >= L]
    random.shuffle(candidates)
    for ps, pe in candidates:
        max_start = pe - L
        if max_start < ps:
            continue
        # try a few random offsets within this piece
        for _ in range(10):
            ws = random.randint(ps, max_start)
            we = ws + L
            if (ws, we) in used:
                continue
            if overlaps_uorf(ws, we):
                continue
            used.add((ws, we))
            negatives.append((ws, we))
            placed = True
            break
        if placed:
            break
    if not placed:
        fails += 1

print(f"Built {len(negatives)} length-matched negatives ({fails} could not be placed)")

with open("uorf_negatives_chr22.bed", "w") as out:
    for ws, we in negatives:
        out.write(f"{CHROM}\t{ws}\t{we}\tneg\t0\t+\n")
print("Wrote uorf_negatives_chr22.bed")