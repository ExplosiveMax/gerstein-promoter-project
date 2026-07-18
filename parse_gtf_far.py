import gzip
import argparse
from gene_overlap import load_gene_intervals, overlaps_any_except

parser = argparse.ArgumentParser()
parser.add_argument("--offset", type=int, default=0)
parser.add_argument("--window", type=int, default=500)
args = parser.parse_args()

GTF_FILE = "gencode.v47.basic.annotation.gtf.gz"
BED_FILE = "promoters.bed"
WINDOW = args.window
HALF = WINDOW // 2
TARGET_COUNT = 500

print("Loading gene intervals for overlap filtering...")
gene_intervals = load_gene_intervals(GTF_FILE, "chr22")
print(f"Loaded {len(gene_intervals)} gene intervals on chr22")

count, skipped, scanned = 0, 0, 0

with gzip.open(GTF_FILE, "rt") as gtf, open(BED_FILE, "w") as bed:
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue
        chrom = fields[0]
        if chrom != "chr22":
            continue
        start, end, strand = int(fields[3]), int(fields[4]), fields[6]
        gene_name = "unknown"
        for part in fields[8].split(";"):
            part = part.strip()
            if part.startswith("gene_name"):
                gene_name = part.split('"')[1]
                break

        scanned += 1
        if strand == "+":
            center = start - args.offset
        else:
            center = end + args.offset
        win_start = max(0, center - HALF)
        win_end = win_start + WINDOW

        if overlaps_any_except(win_start, win_end, gene_name, gene_intervals):
            skipped += 1
            continue

        bed.write(f"{chrom}\t{win_start}\t{win_end}\t{gene_name}\n")
        count += 1
        if count >= TARGET_COUNT:
            break

print(f"Done! Wrote {count} regions ({WINDOW}bp) centered {args.offset}bp upstream of TSS")
print(f"Scanned {scanned} genes, skipped {skipped} for overlapping another gene")
if count < TARGET_COUNT:
    print(f"WARNING: only found {count}/{TARGET_COUNT} — chr22 may be too gene-dense at this offset")