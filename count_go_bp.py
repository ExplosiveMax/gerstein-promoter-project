import gzip
from collections import defaultdict

GAF = "HUMAN-uniprot.gaf.gz"

# gene symbols per GO:BP term
term_genes = defaultdict(set)

with gzip.open(GAF, "rt") as f:
    for line in f:
        if line.startswith("!"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        symbol = c[2]
        go_id = c[4]
        aspect = c[8]
        if aspect != "P":   # Biological Process only
            continue
        term_genes[go_id].add(symbol)

# sort by number of distinct genes
ranked = sorted(term_genes.items(), key=lambda kv: len(kv[1]), reverse=True)
print(f"Total distinct GO:BP terms: {len(ranked)}")
print(f"\nTop 30 BP terms by gene count:")
for go_id, genes in ranked[:30]:
    print(f"  {go_id}\t{len(genes)} genes")
