"""Species registry for the cross-species promoter pipeline.

Everything downstream is keyed off this file, so adding a new species
(chimp, rat, macaque, ...) is a one-entry change, not a code change.

Human is always the reference/query species; every target species is compared
against it via 1-to-1 orthologs pulled from Ensembl BioMart.

IMPORTANT (silent-breakage trap): the BioMart release and the GENCODE
annotation MUST be on the same assembly, or gene IDs will map to the wrong
coordinates without ever erroring.
    human : GRCh38   (GENCODE v47, BioMart current)
    mouse : GRCm39   (GENCODE M39,  BioMart current)   <-- both GRCm39, aligned
    chimp : must confirm current BioMart chimp assembly == your GENCODE/Ensembl FASTA
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    key: str
    label: str                 # 0/1 class label handled at build time; label here is the name
    genome_fasta: str          # soft-masked primary assembly FASTA (has .fai alongside)
    gtf: str                   # matching GENCODE GTF (same assembly!)
    biomart_dataset: str = ""  # only the query species needs this
    homolog_prefix: str = ""   # BioMart homolog attribute prefix for a *target* species


# ---- reference / query species -------------------------------------------------
HUMAN = Species(
    key="human",
    label="human",
    genome_fasta="/Users/maxwellkim/promoter_project/.claude/worktrees/awesome-mendeleev-f66987/GRCh38.primary_assembly.genome.fa",
    gtf="/Users/maxwellkim/promoter_project/.claude/worktrees/awesome-mendeleev-f66987/gencode.v47.basic.annotation.gtf",
    biomart_dataset="hsapiens_gene_ensembl",
)
# ---- target species (compared against human) -----------------------------------
# Fill genome_fasta / gtf after you download GENCODE M39 (GRCm39) for mouse.
# homolog_prefix is the Ensembl species stem used in BioMart homolog attributes,
# e.g. mmusculus -> attribute "mmusculus_homolog_ensembl_gene".
TARGETS = {
    "mouse": Species(
        key="mouse",
        label="mouse",
        genome_fasta="/Users/maxwellkim/promoter_project/mouse/GRCm39.primary_assembly.genome.fa",
        gtf="/Users/maxwellkim/promoter_project/mouse/gencode.vM39.basic.annotation.gtf",
        homolog_prefix="mmusculus",
    ),
    "chimp": Species(
        key="chimp",
        label="chimp",
        genome_fasta="/Users/maxwellkim/promoter_project/chimp/Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa",
        gtf="/Users/maxwellkim/promoter_project/chimp/Pan_troglodytes.Pan_tro_3.0.116.gtf",
        homolog_prefix="ptroglodytes",
    ),
    # Rhesus macaque (~25 My from human) — fills the middle of the divergence
    # curve between chimp (~0.50) and mouse (~0.94). Ensembl release-116 (== the
    # BioMart current release used by get_orthologs.py), assembly Mmul_10.
    "macaque": Species(
        key="macaque",
        label="macaque",
        genome_fasta="/Users/maxwellkim/promoter_project/macaque/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa",
        gtf="/Users/maxwellkim/promoter_project/macaque/Macaca_mulatta.Mmul_10.116.gtf",
        homolog_prefix="mmulatta",
    ),
}