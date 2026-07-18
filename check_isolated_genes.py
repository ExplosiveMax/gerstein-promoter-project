import gzip

GTF = "gencode.v47.basic.annotation.gtf.gz"
MAX_OFFSET = 1000000

genes = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.strip().split("\t")
        if c[2] != "gene" or c[0] != "chr22": continue
        start, end = int(c[3]), int(c[4])
        name = "unknown"
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("gene_name"):
                name = part.split('"')[1]; break
        genes.append((start, end, name))

genes.sort()
isolated = []
for i, (s, e, name) in enumerate(genes):
    ok = True
    for j, (s2, e2, name2) in enumerate(genes):
        if i == j: continue
        if abs(s2 - s) < MAX_OFFSET or abs(e2 - e) < MAX_OFFSET:
            ok = False
            break
    if ok:
        isolated.append((s, e, name))

print(f"Total chr22 genes: {len(genes)}")
print(f"Genes with no neighbor within {MAX_OFFSET}bp: {len(isolated)}")
