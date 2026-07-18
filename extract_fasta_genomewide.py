import argparse
from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument("--genome", default="GRCh38.primary_assembly.genome.fa")
args = parser.parse_args()

genome = Fasta(args.genome)

BEDS = [
    ("functional_promoters_genomewide.bed", "functional_promoters_genomewide.fasta"),
    ("pseudo_promoters_genomewide.bed", "pseudo_promoters_genomewide.fasta"),
    ("random_genomewide.bed", "random_genomewide.fasta"),
]

for bed_path, fasta_path in BEDS:
    n = 0
    with open(bed_path) as bed, open(fasta_path, "w") as out:
        for line in bed:
            fields = line.rstrip("\n").split("\t")
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            seq = str(genome[chrom][s:e]).upper()
            out.write(f">{chrom}:{s}-{e}\n{seq}\n")
            n += 1
    print(f"{bed_path} -> {fasta_path}: {n} sequences")
