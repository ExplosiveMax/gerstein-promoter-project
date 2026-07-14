import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--offset", type=int, default=0)
parser.add_argument("--window", type=int, default=500)
args = parser.parse_args()

# The window is CENTERED at `offset` bp from the TSS.
# offset=0    -> window straddles the TSS
# offset=2000 -> window centered 2000bp upstream of TSS
# Positive offset = upstream (away from gene body), which is where
# we expect signal to decay toward random.

GTF_FILE = "gencode.v47.basic.annotation.gtf.gz"
BED_FILE = "promoters.bed"
WINDOW = args.window
HALF = WINDOW // 2
TEST_LIMIT = 500
count = 0

with gzip.open(GTF_FILE, "rt") as gtf, open(BED_FILE, "w") as bed:
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        if chrom != "chr22":
            continue
        info = fields[8]
        gene_name = "unknown"
        for part in info.split(";"):
            part = part.strip()
            if part.startswith("gene_name"):
                gene_name = part.split('"')[1]
                break
        if strand == "+":
            tss = start
            center = tss - args.offset        # upstream = smaller coord on + strand
        else:
            tss = end
            center = tss + args.offset        # upstream = larger coord on - strand
        win_start = max(0, center - HALF)
        win_end = win_start + WINDOW
        bed.write(f"{chrom}\t{win_start}\t{win_end}\t{gene_name}\n")
        count += 1
        if count >= TEST_LIMIT:
            break

print(f"Done! Wrote {count} regions ({WINDOW}bp) centered {args.offset}bp upstream of TSS")