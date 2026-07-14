import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--window", type=int, default=1000)
args = parser.parse_args()

GTF_FILE = "gencode.v47.basic.annotation.gtf.gz"
BED_FILE = "promoters.bed"
UPSTREAM = args.window
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
            prom_start = max(0, start - UPSTREAM)
            prom_end = start
        else:
            prom_start = end
            prom_end = end + UPSTREAM
        bed.write(f"{chrom}\t{prom_start}\t{prom_end}\t{gene_name}\n")
        count += 1
        if count >= TEST_LIMIT:
            break

print(f"Done! Wrote {count} promoter regions to {BED_FILE} with window size {UPSTREAM}")