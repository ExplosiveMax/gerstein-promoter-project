import gzip
from collections import defaultdict

GAF = "HUMAN-uniprot.gaf.gz"
CATEGORIES = {
    "GO:0006955": "immune_response",
    "GO:0007186": "gpcr_signaling",
    "GO:0055085": "transmembrane_transport",
}

gene_cats = defaultdict(set)
with gzip.open(GAF, "rt") as f:
    for line in f:
        if line.startswith("!"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        symbol, go_id, aspect = c[2], c[4], c[8]
        if aspect != "P": continue
        if go_id in CATEGORIES:
            gene_cats[symbol].add(go_id)

clean = defaultdict(list)
multi = 0
for symbol, cats in gene_cats.items():
    if len(cats) == 1:
        clean[CATEGORIES[next(iter(cats))]].append(symbol)
    else:
        multi += 1

print(f"Multi-category excluded: {multi}")
for name in CATEGORIES.values():
    print(f"  {name}: {len(clean[name])}")
for name, genes in clean.items():
    with open(f"go_{name}_genes.txt", "w") as out:
        for g in sorted(genes):
            out.write(g + "\n")
