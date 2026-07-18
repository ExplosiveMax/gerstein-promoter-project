import random
import argparse

from pyfaidx import Fasta
from gene_overlap_genomewide import (
    STANDARD_CHROMS,
    load_gene_intervals_genomewide,
    overlaps_any,
)

parser = argparse.ArgumentParser()
parser.add_argument("--genome", default="GRCh38.primary_assembly.genome.fa")
parser.add_argument("--gc_tol", type=float, default=0.10)  # matches gcmatch_window_v3.py
parser.add_argument("--max_repeat", type=float, default=0.25)  # matches the rest of the pipeline
parser.add_argument("--max_n", type=float, default=0.05)  # reject assembly-gap-heavy windows
parser.add_argument("--max_attempts", type=int, default=2000)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
random.seed(args.seed)


def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    total = gc + at
    return gc / total if total else 0


def repeat_frac(seq):
    return sum(1 for c in seq if c.islower()) / len(seq)


def n_frac(seq):
    return seq.upper().count("N") / len(seq)


print(f"Opening indexed genome {args.genome} (pyfaidx, no full load into RAM)...")
genome = Fasta(args.genome)
chrom_lens = {c: len(genome[c]) for c in STANDARD_CHROMS if c in genome.keys()}
print(f"Indexed {len(chrom_lens)} standard chromosomes")

print("Loading protein-coding + pseudogene intervals for overlap exclusion...")
protein_iv = load_gene_intervals_genomewide(
    "gencode.v47.basic.annotation.gtf.gz", "gene_intervals_genomewide_protein_coding.pkl", feature="gene"
)
pseudo_iv = load_gene_intervals_genomewide(
    "gencode.v47.2wayconspseudos.gtf.gz", "gene_intervals_genomewide_2way_pseudo.pkl", feature="transcript"
)
combined_iv = {c: sorted(protein_iv.get(c, []) + pseudo_iv.get(c, [])) for c in STANDARD_CHROMS}

positives = []
for bedfile in ("functional_promoters_genomewide.bed", "pseudo_promoters_genomewide.bed"):
    with open(bedfile) as f:
        for line in f:
            chrom, s, e, name = line.rstrip("\n").split("\t")
            positives.append((chrom, int(s), int(e)))
print(f"Matching {len(positives)} positives (functional + pseudogene combined)")

chroms_list = list(chrom_lens.keys())
regions = []
fails = 0
for chrom, s, e in positives:
    window = e - s
    target_gc = gc_content(str(genome[chrom][s:e]))
    placed = False
    for _ in range(args.max_attempts):
        rchrom = random.choice(chroms_list)
        clen = chrom_lens[rchrom]
        if clen <= window:
            continue
        start = random.randint(0, clen - window)
        end = start + window
        if overlaps_any(start, end, combined_iv[rchrom]):
            continue
        seq = str(genome[rchrom][start:end])
        if repeat_frac(seq) > args.max_repeat:
            continue
        if n_frac(seq) > args.max_n:
            continue
        if abs(gc_content(seq) - target_gc) <= args.gc_tol:
            regions.append((rchrom, start, end))
            placed = True
            break
    if not placed:
        fails += 1

print(f"Found {len(regions)} GC-matched (per-sequence, +/-{args.gc_tol}), non-genic regions genome-wide")
print(f"Failed to place: {fails}")

with open("random_genomewide.bed", "w") as f:
    for chrom, start, end in regions:
        f.write(f"{chrom}\t{start}\t{end}\n")
print("Saved random_genomewide.bed")
