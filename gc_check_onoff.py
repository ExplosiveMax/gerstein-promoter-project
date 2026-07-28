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
    g = seq.count("G")+seq.count("C"); at = seq.count("A")+seq.count("T")
    return g/(g+at) if (g+at) else 0

for tissue in ["Liver", "Muscle_Skeletal", "Whole_Blood"]:
    for label in ["high", "low"]:
        s = load_fasta(f"onoff_{tissue}_{label}.fasta")
        g = [gc(x) for x in s]
        print(f"{tissue} {label}: mean GC {sum(g)/len(g):.4f} (n={len(s)})")
