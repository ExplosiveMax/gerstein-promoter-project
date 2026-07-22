import argparse
from pyfaidx import Fasta
parser = argparse.ArgumentParser()
parser.add_argument("--genome", default="GRCh38.primary_assembly.genome.fa")
args = parser.parse_args()
genome = Fasta(args.genome)
BEDS = [
    ("go_immune_response_promoters.bed", "go_immune_response.fasta"),
    ("go_gpcr_signaling_promoters.bed", "go_gpcr_signaling.fasta"),
    ("go_transmembrane_transport_promoters.bed", "go_transmembrane_transport.fasta"),
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
