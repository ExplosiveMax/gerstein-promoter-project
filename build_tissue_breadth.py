import gzip

GTEX = "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
TPM_THRESHOLD = 1.0      # gene counts as "expressed" in a tissue if median TPM >= this
BROAD_MIN = 64           # expressed in >= this many tissues -> housekeeping/broad
SPECIFIC_MAX = 5         # expressed in <= this many (but >=1) -> tissue-specific

housekeeping = []
tissue_specific = []
n_tissues = None

with gzip.open(GTEX, "rt") as f:
    f.readline()                      # skip "#1.2"
    dims = f.readline().split()       # "74628  68"
    header = f.readline().rstrip("\n").split("\t")
    n_tissues = len(header) - 2       # first two cols are Name, Description
    print(f"Tissues: {n_tissues}, genes reported: {dims[0]}")

    for line in f:
        parts = line.rstrip("\n").split("\t")
        gene_id = parts[0]
        tpms = parts[2:]
        n_expressed = sum(1 for v in tpms if float(v) >= TPM_THRESHOLD)
        if n_expressed >= BROAD_MIN:
            housekeeping.append((gene_id, n_expressed))
        elif 1 <= n_expressed <= SPECIFIC_MAX:
            tissue_specific.append((gene_id, n_expressed))

print(f"Housekeeping (>={BROAD_MIN}/{n_tissues} tissues): {len(housekeeping)}")
print(f"Tissue-specific (1-{SPECIFIC_MAX} tissues):        {len(tissue_specific)}")

with open("housekeeping_genes.txt", "w") as out:
    for gid, n in housekeeping:
        out.write(f"{gid}\t{n}\n")
with open("tissuespecific_genes.txt", "w") as out:
    for gid, n in tissue_specific:
        out.write(f"{gid}\t{n}\n")
print("Wrote housekeeping_genes.txt and tissuespecific_genes.txt")
