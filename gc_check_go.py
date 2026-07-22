def load_fasta(path):
    seqs, cur = [], ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else: cur += line
    if cur: seqs.append(cur.upper())
    return seqs

def gc(seq):
    g = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return g/(g+at) if (g+at) else 0

for name, path in [("immune","go_immune_response.fasta"),
                   ("gpcr","go_gpcr_signaling.fasta"),
                   ("transmembrane","go_transmembrane_transport.fasta")]:
    s = load_fasta(path)
    g = [gc(x) for x in s]
    print(f"{name}: mean GC {sum(g)/len(g):.4f}  (n={len(s)})")
