import random
random.seed(42)

def load_fasta(path):
    seqs, cur, hdr, out = [], "", None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur: out.append((hdr, cur.upper()))
            hdr, cur = line, ""
        else:
            cur += line
    if cur: out.append((hdr, cur.upper()))
    return out

def gc(seq):
    g = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return g / (g + at) if (g + at) else 0

hk = load_fasta("housekeeping_promoters.fasta")
ts = load_fasta("tissuespecific_promoters.fasta")

BIN = 0.02  # GC bin width
def binkey(seq): return round(gc(seq) / BIN)

from collections import defaultdict
hk_bins, ts_bins = defaultdict(list), defaultdict(list)
for h, s in hk: hk_bins[binkey(s)].append((h, s))
for h, s in ts: ts_bins[binkey(s)].append((h, s))

# For each GC bin, keep min(count_hk, count_ts) from each class
hk_matched, ts_matched = [], []
for b in sorted(set(hk_bins) | set(ts_bins)):
    n = min(len(hk_bins[b]), len(ts_bins[b]))
    if n == 0: continue
    random.shuffle(hk_bins[b]); random.shuffle(ts_bins[b])
    hk_matched += hk_bins[b][:n]
    ts_matched += ts_bins[b][:n]

def write(path, items):
    with open(path, "w") as f:
        for h, s in items: f.write(f"{h}\n{s}\n")

write("housekeeping_gcmatched.fasta", hk_matched)
write("tissuespecific_gcmatched.fasta", ts_matched)

hk_gc = sum(gc(s) for _, s in hk_matched) / len(hk_matched)
ts_gc = sum(gc(s) for _, s in ts_matched) / len(ts_matched)
print(f"Matched: housekeeping {len(hk_matched)}, tissue-specific {len(ts_matched)}")
print(f"Mean GC after matching: housekeeping {hk_gc:.4f}, tissue-specific {ts_gc:.4f}")
print(f"(before: 0.520 vs 0.462; want these two now nearly equal)")
