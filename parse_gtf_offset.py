import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--offset", type=int, default=0)
args = parser.parse_args()

# offset=0: window is [-1000, 0] relative to TSS (current default)
# offset=+500: window is [-500, +500]
# offset=+1000: window is [0, +1000] (into gene body)
# offset=-500: window is [-1500, -500] (further upstream)
# offset=-1000: window is [-2000, -1000]

GTF_FILE = "gencode.v47.basic.annotation.gtf.gz"
BED_FILE = "promoters.bed"
WINDOW = 1000
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
            win_end = tss + args.offset
            win_start = max(0, win_end - WINDOW)
        else:
            tss = end
            win_start = tss - args.offset
            win_end = win_start + WINDOW
        bed.write(f"{chrom}\t{win_start}\t{win_end}\t{gene_name}\n")
        count += 1
        if count >= TEST_LIMIT:
            break

print(f"Done! Wrote {count} regions to {BED_FILE} with offset {args.offset}")