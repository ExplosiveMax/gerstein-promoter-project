def load_fasta(path):
    seqs, cur = [], ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs

def gc(seq):
    g = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return g / (g + at) if (g + at) else 0

def cpg_oe(seq):
    # observed/expected CpG ratio: the standard CpG-island metric
    c, g = seq.count("C"), seq.count("G")
    cg = seq.count("CG")
    L = len(seq)
    if c == 0 or g == 0: return 0
    expected = (c * g) / L
    return cg / expected if expected else 0

for name, path in [("housekeeping", "housekeeping_promoters.fasta"),
                   ("tissue-specific", "tissuespecific_promoters.fasta")]:
    seqs = load_fasta(path)
    gcs = [gc(s) for s in seqs]
    oes = [cpg_oe(s) for s in seqs]
    print(f"{name} (n={len(seqs)}):")
    print(f"  mean GC:      {sum(gcs)/len(gcs):.4f}")
    print(f"  mean CpG o/e: {sum(oes)/len(oes):.4f}")
    # fraction that qualify as CpG-island-like (GC>0.5 and o/e>0.6, standard heuristic)
    island = sum(1 for g,o in zip(gcs,oes) if g>0.5 and o>0.6)
    print(f"  CpG-island-like: {island}/{len(seqs)} ({100*island/len(seqs):.0f}%)")
