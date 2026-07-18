import gzip

GTF = "gencode.v47.basic.annotation.gtf.gz"
MAX_OFFSET = 10000

genes = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.strip().split("\t")
        if c[2] != "gene" or c[0] != "chr22": continue
        start, end, strand = int(c[3]), int(c[4]), c[6]
        name = "unknown"
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("gene_name"):
                name = part.split('"')[1]; break
        genes.append((start, end, strand, name))

genes.sort()
isolated = []
for i, (s, e, strand, name) in enumerate(genes):
    ok = True
    for j, (s2, e2, strand2, name2) in enumerate(genes):
        if i == j: continue
        if abs(s2 - s) < MAX_OFFSET or abs(e2 - e) < MAX_OFFSET:
            ok = False; break
    if ok:
        isolated.append((s, e, strand, name))

with open("isolated_genes_10kb.tsv", "w") as f:
    for s, e, strand, name in isolated:
        f.write(f"{s}\t{e}\t{strand}\t{name}\n")
print(f"Saved {len(isolated)} isolated genes to isolated_genes_10kb.tsv")
