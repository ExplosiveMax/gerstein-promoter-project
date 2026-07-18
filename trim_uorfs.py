def load_fasta(filepath):
    seqs, cur = [], ""
    for line in open(filepath):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs

seqs = load_fasta("uorfs_capped.fasta")
trimmed = [s[3:-3] for s in seqs if len(s) > 12]
print(f"{len(seqs)} -> {len(trimmed)} (dropped {len(seqs)-len(trimmed)} too short to trim)")

with open("uorfs_capped_trimmed.fasta", "w") as out:
    for i, s in enumerate(trimmed):
        out.write(f">trimmed_{i}\n{s}\n")
