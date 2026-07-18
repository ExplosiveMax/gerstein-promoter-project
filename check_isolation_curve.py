import gzip

GTF = "gencode.v47.basic.annotation.gtf.gz"

genes = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.strip().split("\t")
        if c[2] != "gene" or c[0] != "chr22": continue
        start, end = int(c[3]), int(c[4])
        genes.append((start, end))

genes.sort()

for max_offset in [10000, 25000, 50000, 100000, 150000, 200000, 300000, 500000]:
    isolated = 0
    for i, (s, e) in enumerate(genes):
        ok = True
        for j, (s2, e2) in enumerate(genes):
            if i == j: continue
            if abs(s2 - s) < max_offset or abs(e2 - e) < max_offset:
                ok = False
                break
        if ok:
            isolated += 1
    print(f"max_offset={max_offset}: {isolated} isolated genes")
