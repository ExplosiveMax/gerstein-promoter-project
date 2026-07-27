import gzip
import random
import argparse
from gene_overlap_genomewide import (
    STANDARD_CHROMS,
    load_gene_intervals_genomewide,
    overlaps_any_except,
)

parser = argparse.ArgumentParser()
parser.add_argument("--tissue", required=True)
parser.add_argument("--window", type=int, default=1000)
parser.add_argument("--limit", type=int, default=500)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
random.seed(args.seed)

BASIC_GTF = "gencode.v47.basic.annotation.gtf.gz"

def parse_attr(info, key):
    for part in info.split(";"):
        part = part.strip()
        if part.startswith(key):
            return part.split('"')[1]
    return "unknown"

def promoter_coords(start, end, strand, window):
    if strand == "+":
        return max(0, start - window), start
    else:
        return end, end + window

def load_ids(path):
    return set(l.strip() for l in open(path) if l.strip())

high_ids = load_ids(f"onoff_{args.tissue}_high.txt")
low_ids  = load_ids(f"onoff_{args.tissue}_low.txt")
print(f"{args.tissue}: high={len(high_ids)}, low={len(low_ids)}")

print("Loading gene intervals per chromosome (for overlap exclusion)...")
gene_intervals = load_gene_intervals_genomewide(
    BASIC_GTF, "gene_intervals_genomewide_protein_coding.pkl", feature="gene"
)

found = {}
all_wanted = high_ids | low_ids
with gzip.open(BASIC_GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene" or c[0] not in STANDARD_CHROMS: continue
        gid = parse_attr(c[8], "gene_id")
        if gid in all_wanted and gid not in found:
            found[gid] = (c[0], int(c[3]), int(c[4]), c[6], parse_attr(c[8], "gene_name"))

print(f"Located in GTF: {len(found)}/{len(all_wanted)}")

def build(ids, label, out_path):
    cands = [gid for gid in ids if gid in found]
    random.shuffle(cands)
    written, skipped = 0, 0
    with open(out_path, "w") as out:
        for gid in cands:
            if written >= args.limit: break
            chrom, start, end, strand, gname = found[gid]
            ps, pe = promoter_coords(start, end, strand, args.window)
            chrom_intervals = gene_intervals.get(chrom, [])
            if overlaps_any_except(ps, pe, gname, chrom_intervals):
                skipped += 1; continue
            out.write(f"{chrom}\t{ps}\t{pe}\t{gid}\t0\t{strand}\n")
            written += 1
    print(f"{label}: wrote {written} promoters ({skipped} skipped overlap, {len(cands)} candidates)")

build(high_ids, f"{args.tissue}_high", f"onoff_{args.tissue}_high.bed")
build(low_ids, f"{args.tissue}_low", f"onoff_{args.tissue}_low.bed")
