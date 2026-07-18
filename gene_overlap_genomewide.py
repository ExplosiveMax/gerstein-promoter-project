import gzip
import os
import pickle

STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def _parse_attr(info, key):
    for part in info.split(";"):
        part = part.strip()
        if part.startswith(key):
            return part.split('"')[1]
    return "unknown"


def load_gene_intervals_genomewide(gtf_file, cache_file, feature="gene", chroms=STANDARD_CHROMS):
    """Load (start, end, name) intervals per chromosome from a GTF/GFF-style file.

    Caches the parsed result as a pickle next to cache_file so repeated runs
    (e.g. across offset/window sweeps) don't have to re-parse a multi-GB gtf.
    """
    if os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            return pickle.load(f)

    intervals = {c: [] for c in chroms}
    opener = gzip.open if gtf_file.endswith(".gz") else open
    with opener(gtf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != feature or fields[0] not in chroms:
                continue
            start, end = int(fields[3]), int(fields[4])
            name = _parse_attr(fields[8], "gene_name")
            intervals[fields[0]].append((start, end, name))

    for c in intervals:
        intervals[c].sort()

    with open(cache_file, "wb") as f:
        pickle.dump(intervals, f)
    return intervals


def overlaps_any(win_start, win_end, intervals_for_chrom):
    for s, e, _ in intervals_for_chrom:
        if s >= win_end:
            break
        if e > win_start:
            return True
    return False


def overlaps_any_except(win_start, win_end, self_name, intervals_for_chrom):
    for s, e, name in intervals_for_chrom:
        if s >= win_end:
            break
        if e > win_start and name != self_name:
            return True
    return False
