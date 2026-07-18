import gzip
import random
import argparse
from collections import defaultdict

from gene_overlap_genomewide import (
    STANDARD_CHROMS,
    load_gene_intervals_genomewide,
    overlaps_any_except,
)

parser = argparse.ArgumentParser()
parser.add_argument("--window", type=int, default=1000)
parser.add_argument("--limit", type=int, default=300)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
random.seed(args.seed)

BASIC_GTF = "gencode.v47.basic.annotation.gtf.gz"
PSEUDO_GTF = "gencode.v47.2wayconspseudos.gtf.gz"
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


print("Loading protein-coding gene intervals genome-wide (for overlap exclusion + functional promoters)...")
gene_intervals = load_gene_intervals_genomewide(
    BASIC_GTF, "gene_intervals_genomewide_protein_coding.pkl", feature="gene"
)
total_genes = sum(len(v) for v in gene_intervals.values())
print(f"Loaded {total_genes} gene intervals across {len(gene_intervals)} chromosomes")

# Functional candidates: protein_coding genes, re-parsed to keep gene_type + strand
functional_candidates = []
with gzip.open(BASIC_GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene" or c[0] not in STANDARD_CHROMS:
            continue
        gtype = parse_attr(c[8], "gene_type")
        if gtype != "protein_coding":
            continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        name = parse_attr(c[8], "gene_name")
        functional_candidates.append((c[0], start, end, strand, name))

print(f"Protein-coding gene candidates genome-wide: {len(functional_candidates)}")
random.shuffle(functional_candidates)

# Pseudogene candidates from the Yale-UCSC 2-way consensus set (transcript records only)
pseudo_candidates = []
with gzip.open(PSEUDO_GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "transcript" or c[0] not in STANDARD_CHROMS:
            continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        name = parse_attr(c[8], "gene_id")
        pseudo_candidates.append((c[0], start, end, strand, name))

print(f"Yale-UCSC 2-way consensus pseudogene candidates genome-wide: {len(pseudo_candidates)}")
random.shuffle(pseudo_candidates)


def build(candidates, limit, label):
    kept = []
    skipped_overlap = 0
    skipped_boundary = 0
    for chrom, start, end, strand, name in candidates:
        if len(kept) >= limit:
            break
        ps, pe = promoter_coords(start, end, strand, WINDOW)
        if ps < 0:
            skipped_boundary += 1
            continue
        # exclude promoter windows that fall inside another annotated gene body
        # (self-overlap is fine for functional genes; pseudogene names don't
        # collide with gene_name space so "except self" is a no-op there)
        if overlaps_any_except(ps, pe, name, gene_intervals[chrom]):
            skipped_overlap += 1
            continue
        kept.append((chrom, ps, pe, name))
    print(f"{label}: kept {len(kept)} (skipped {skipped_overlap} gene-overlap, {skipped_boundary} boundary)")
    return kept


functional = build(functional_candidates, LIMIT, "functional")
pseudo = build(pseudo_candidates, LIMIT, "pseudogene")

with open("functional_promoters_genomewide.bed", "w") as out:
    for chrom, s, e, name in functional:
        out.write(f"{chrom}\t{s}\t{e}\t{name}\n")
with open("pseudo_promoters_genomewide.bed", "w") as out:
    for chrom, s, e, name in pseudo:
        out.write(f"{chrom}\t{s}\t{e}\t{name}\n")

print("Wrote functional_promoters_genomewide.bed and pseudo_promoters_genomewide.bed")
