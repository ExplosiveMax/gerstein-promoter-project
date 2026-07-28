from pyfaidx import Fasta

genome = Fasta("/Users/maxwellkim/promoter_project/.claude/worktrees/awesome-mendeleev-f66987/GRCh38.primary_assembly.genome.fa")
CATEGORIES = ["signal_transduction", "protein_ubiquitination", "proteolysis"]

for name in CATEGORIES:
    bed_path = f"go_{name}_promoters.bed"
    fasta_path = f"go_{name}_pw.fasta"
    n = 0
    with open(bed_path) as bed, open(fasta_path, "w") as out:
        for line in bed:
            fields = line.rstrip("\n").split("\t")
            chrom, s, e = fields[0], int(fields[1]), int(fields[2])
            seq = str(genome[chrom][s:e]).upper()
            out.write(f">{chrom}:{s}-{e}\n{seq}\n")
            n += 1
    print(f"{bed_path} -> {fasta_path}: {n} sequences")
