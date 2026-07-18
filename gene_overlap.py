import gzip
import os

def load_gene_intervals(gtf_file, chrom, cache_file="chr22_gene_intervals.tsv"):
    if os.path.exists(cache_file):
        intervals = []
        with open(cache_file) as f:
            for line in f:
                s, e, name = line.rstrip("\n").split("\t")
                intervals.append((int(s), int(e), name))
        return intervals

    intervals = []
    with gzip.open(gtf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "gene" or fields[0] != chrom:
                continue
            start, end = int(fields[3]), int(fields[4])
            name = "unknown"
            for part in fields[8].split(";"):
                part = part.strip()
                if part.startswith("gene_name"):
                    name = part.split('"')[1]
                    break
            intervals.append((start, end, name))
    intervals.sort()

    with open(cache_file, "w") as f:
        for s, e, name in intervals:
            f.write(f"{s}\t{e}\t{name}\n")
    return intervals

def overlaps_any_except(win_start, win_end, self_name, intervals):
    for s, e, name in intervals:
        if s >= win_end:
            break
        if e > win_start and name != self_name:
            return True
    return False

def overlaps_any(win_start, win_end, intervals):
    for s, e, name in intervals:
        if s >= win_end:
            break
        if e > win_start:
            return True
    return False