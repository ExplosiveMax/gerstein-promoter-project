import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--offset", type=int, default=0)
parser.add_argument("--window", type=int, default=500)
args = parser.parse_args()

CHR_LEN = 50818468
HALF = args.window // 2
count = 0

with open("isolated_genes_10kb.tsv") as f, open("promoters.bed", "w") as bed:
    for line in f:
        s, e, strand, name = line.rstrip("\n").split("\t")
        s, e = int(s), int(e)
        if strand == "+":
            center = s - args.offset
        else:
            center = e + args.offset
        win_start = max(0, center - HALF)
        win_end = win_start + args.window
        if win_end > CHR_LEN or win_start < 0:
            continue
        bed.write(f"chr22\t{win_start}\t{win_end}\t{name}\n")
        count += 1

print(f"Wrote {count} regions at offset {args.offset}")
