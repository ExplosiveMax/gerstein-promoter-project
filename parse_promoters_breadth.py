import gzip
import random
import argparse
from gene_overlap_genomewide import (
    STANDARD_CHROMS,
    load_gene_intervals_genomewide,
    overlaps_any_except,
)

parser = argparse.ArgumentParser()
parser.add_argument("--window", type=int, default=1000)
parser.add_argument("--limit", type=int, default=500)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
random.seed(args.seed)

BASIC_GTF = "gencode.v47.basic.annotation.gtf.gz"
WINDOW = args.window
LIMIT = args.limit

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

def load_id_set(path):
    ids = set()
    with open(path) as f:
        for line in f:
            ids.add(line.split("\t")[0].strip())
    return ids

print("Loading gene-ID class lists...")
housekeeping_ids = load_id_set("housekeeping_genes.txt")
specific_ids     = load_id_set("tissuespecific_genes.txt")
print(f"Housekeeping IDs: {len(housekeeping_ids)}, tissue-specific IDs: {len(specific_ids)}")

print("Loading gene intervals per chromosome (for overlap exclusion)...")
gene_intervals = load_gene_intervals_genomewide(
    BASIC_GTF, "gene_intervals_genomewide_protein_coding.pkl", feature="gene"
)

def collect(gtf, want_ids):
    found = {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "gene" or c[0] not in STANDARD_CHROMS: continue
            gid = parse_attr(c[8], "gene_id")
            if gid in want_ids and gid not in found:
                found[gid] = (c[0], int(c[3]), int(c[4]), c[6], parse_attr(c[8], "gene_name"))
    return found

print("Scanning GTF for class members...")
hk_genes = collect(BASIC_GTF, housekeeping_ids)
ts_genes = collect(BASIC_GTF, specific_ids)
print(f"Located in GTF: housekeeping {len(hk_genes)}, tissue-specific {len(ts_genes)}")

def build_promoters(genes, label, out_path):
    items = list(genes.items())
    random.shuffle(items)
    written = 0
    skipped_overlap = 0
    with open(out_path, "w") as out:
        for gid, (chrom, start, end, strand, gname) in items:
            if written >= LIMIT: break
            ps, pe = promoter_coords(start, end, strand, WINDOW)
            chrom_intervals = gene_intervals.get(chrom, [])
            if overlaps_any_except(ps, pe, gname, chrom_intervals):
                skipped_overlap += 1
                continue
            out.write(f"{chrom}\t{ps}\t{pe}\t{gid}\t0\t{strand}\n")
            written += 1
    print(f"{label}: wrote {written} promoters ({skipped_overlap} skipped for gene overlap)")

build_promoters(hk_genes, "housekeeping", "housekeeping_promoters.bed")
build_promoters(ts_genes, "tissuespecific", "tissuespecific_promoters.bed")
