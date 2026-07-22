import random
from collections import defaultdict
random.seed(42)

def load_fasta(path):
    out, cur, hdr = [], "", None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: out.append((hdr, cur.upper()))
            hdr, cur = line, ""
        else: cur += line
    if cur: out.append((hdr, cur.upper()))
    return out

def gc(seq):
    g = seq.count("G")+seq.count("C"); at = seq.count("A")+seq.count("T")
    return g/(g+at) if (g+at) else 0

BIN = 0.02
classes = {
    "immune": load_fasta("go_immune_response.fasta"),
    "gpcr": load_fasta("go_gpcr_signaling.fasta"),
    "transmembrane": load_fasta("go_transmembrane_transport.fasta"),
}

binned = {name: defaultdict(list) for name in classes}
for name, seqs in classes.items():
    for h, s in seqs:
        binned[name][round(gc(s)/BIN)].append((h, s))

all_bins = set()
for name in classes: all_bins |= set(binned[name])

matched = {name: [] for name in classes}
for b in sorted(all_bins):
    counts = [len(binned[name][b]) for name in classes]
    n = min(counts)
    if n == 0: continue
    for name in classes:
        random.shuffle(binned[name][b])
        matched[name] += binned[name][b][:n]

outnames = {"immune":"go_immune_gcmatched.fasta",
            "gpcr":"go_gpcr_gcmatched.fasta",
            "transmembrane":"go_transmembrane_gcmatched.fasta"}
for name, items in matched.items():
    with open(outnames[name], "w") as f:
        for h, s in items: f.write(f"{h}\n{s}\n")
    g = sum(gc(s) for _,s in items)/len(items)
    print(f"{name}: {len(items)} kept, mean GC {g:.4f}")
