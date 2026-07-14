import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--window", type=int, default=1000)
parser.add_argument("--limit", type=int, default=300)
args = parser.parse_args()

GTF = "gencode.v47.basic.annotation.gtf.gz"
CHROM = "chr22"
WINDOW = args.window
LIMIT = args.limit

PSEUDO_TYPES = {
    "processed_pseudogene", "unprocessed_pseudogene",
    "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
    "transcribed_unitary_pseudogene", "unitary_pseudogene",
    "rRNA_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene",
    "TR_V_pseudogene", "TR_J_pseudogene", "IG_J_pseudogene",
    "translated_processed_pseudogene", "IG_pseudogene",
}

def gene_type_of(info):
    for part in info.split(";"):
        part = part.strip()
        if part.startswith("gene_type"):
            return part.split('"')[1]
    return None

def promoter_coords(start, end, strand, window):
    if strand == "+":
        return max(0, start - window), start
    else:
        return end, end + window

functional = []
pseudo = []

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] != CHROM or c[2] != "gene":
            continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        gtype = gene_type_of(c[8])
        ps, pe = promoter_coords(start, end, strand, WINDOW)
        if gtype == "protein_coding":
            functional.append((ps, pe))
        elif gtype in PSEUDO_TYPES:
            pseudo.append((ps, pe))

print(f"Functional promoters available: {len(functional)}")
print(f"Pseudogene promoters available: {len(pseudo)}")

functional = functional[:LIMIT]
pseudo = pseudo[:LIMIT]
print(f"Using {len(functional)} functional, {len(pseudo)} pseudogene")

with open("functional_promoters.bed", "w") as out:
    for s, e in functional:
        out.write(f"{CHROM}\t{s}\t{e}\n")
with open("pseudo_promoters.bed", "w") as out:
    for s, e in pseudo:
        out.write(f"{CHROM}\t{s}\t{e}\n")

print("Wrote functional_promoters.bed and pseudo_promoters.bed")