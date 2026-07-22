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
parser.add_argument("--limit", type=int, default=300)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
random.seed(args.seed)

BASIC_GTF = "gencode.v47.basic.annotation.gtf.gz"

CATEGORIES = ["immune_response", "gpcr_signaling", "transmembrane_transport"]

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

def load_symbols(name):
    with open(f"go_{name}_genes.txt") as f:
        return set(l.strip() for l in f if l.strip())

# symbol sets per category
cat_symbols = {name: load_symbols(name) for name in CATEGORIES}
all_wanted = set().union(*cat_symbols.values())
print("Loaded symbol lists:", {k: len(v) for k, v in cat_symbols.items()})

print("Loading gene intervals per chromosome (for overlap exclusion)...")
gene_intervals = load_gene_intervals_genomewide(
    BASIC_GTF, "gene_intervals_genomewide_protein_coding.pkl", feature="gene"
)

# scan GTF: map gene_name (symbol) -> coords, keep only wanted symbols
found = {}  # symbol -> (chrom,start,end,strand,gene_id)
with gzip.open(BASIC_GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene" or c[0] not in STANDARD_CHROMS: continue
        gname = parse_attr(c[8], "gene_name")
        if gname in all_wanted and gname not in found:
            found[gname] = (c[0], int(c[3]), int(c[4]), c[6], parse_attr(c[8], "gene_id"))

print(f"Symbols located in GTF: {len(found)}/{len(all_wanted)}")

def build(name, out_path, label):
    syms = [s for s in cat_symbols[name] if s in found]
    random.shuffle(syms)
    written, skipped = 0, 0
    with open(out_path, "w") as out:
        for s in syms:
            if written >= args.limit: break
            chrom, start, end, strand, gid = found[s]
            ps, pe = promoter_coords(start, end, strand, args.window)
            chrom_intervals = gene_intervals.get(chrom, [])
            if overlaps_any_except(ps, pe, s, chrom_intervals):
                skipped += 1
                continue
            out.write(f"{chrom}\t{ps}\t{pe}\t{gid}|{name}\t0\t{strand}\n")
            written += 1
    print(f"{name}: wrote {written} promoters ({skipped} skipped overlap, {len(syms)} candidates)")

for i, name in enumerate(CATEGORIES):
    build(name, f"go_{name}_promoters.bed", i)
