import gzip

GTEX = "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
TISSUES = ["Liver", "Muscle_Skeletal", "Whole_Blood"]
N_PER_GROUP = 1500

with gzip.open(GTEX, "rt") as f:
    f.readline()  # #1.2
    f.readline()  # dims
    header = f.readline().rstrip("\n").split("\t")
    tissue_idx = {t: header.index(t) for t in TISSUES}
    print("Column indices:", tissue_idx)

    rows = []  # (gene_id, [tpm per target tissue])
    for line in f:
        parts = line.rstrip("\n").split("\t")
        gene_id = parts[0]
        vals = {t: float(parts[tissue_idx[t]]) for t in TISSUES}
        rows.append((gene_id, vals))

print(f"Total genes: {len(rows)}")

for t in TISSUES:
    ranked = sorted(rows, key=lambda r: r[1][t], reverse=True)
    top = ranked[:N_PER_GROUP]
    bottom = ranked[-N_PER_GROUP:]
    print(f"{t}: top TPM range {top[-1][1][t]:.2f}-{top[0][1][t]:.2f}, "
          f"bottom TPM range {bottom[0][1][t]:.2f}-{bottom[-1][1][t]:.2f}")
    with open(f"onoff_{t}_high.txt", "w") as out:
        for gid, _ in top: out.write(gid + "\n")
    with open(f"onoff_{t}_low.txt", "w") as out:
        for gid, _ in bottom: out.write(gid + "\n")
print("Wrote onoff_<tissue>_high.txt / _low.txt for each tissue")
