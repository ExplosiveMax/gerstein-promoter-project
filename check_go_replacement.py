import gzip
from collections import defaultdict

GAF = "HUMAN-uniprot.gaf.gz"
CANDIDATES = {
    "GO:0007165": "signal_transduction",
    "GO:0016567": "protein_ubiquitination",
    "GO:0006508": "proteolysis",
}

gene_cats = defaultdict(set)
with gzip.open(GAF, "rt") as f:
    for line in f:
        if line.startswith("!"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        symbol, go_id, aspect = c[2], c[4], c[8]
        if aspect != "P": continue
        if go_id in CANDIDATES:
            gene_cats[symbol].add(go_id)

clean = defaultdict(list)
for symbol, cats in gene_cats.items():
    if len(cats) == 1:
        clean[CANDIDATES[next(iter(cats))]].append(symbol)

# quick match check against GTF gene_names
import gzip as gz
gtf_names = set()
with gz.open("gencode.v47.basic.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene": continue
        for part in c[8].split(";"):
            part = part.strip()
            if part.startswith("gene_name"):
                gtf_names.add(part.split('"')[1])

for name, syms in clean.items():
    matched = sum(1 for s in syms if s in gtf_names)
    print(f"{name}: {len(syms)} genes, {matched} match GTF ({100*matched/len(syms):.0f}%)")
