"""Build a human-vs-target cross-species promoter classification dataset.

Design goals (this is the whole reason the experiment is defensible):

  1. MATCHED GENE SET. For every human promoter we use the promoter of that
     gene's 1:1 ortholog in the target species. Same biological functions on
     both sides -> any separability is species sequence divergence, not
     "human happened to sample immune genes".

  2. REPEAT CONTROL (the big cross-species trap). Human and mouse carry
     DIFFERENT transposon families (Alu vs B1/B2). A classifier can trivially
     "detect species" by spotting lineage-specific repeats instead of promoter
     biology. We use the soft-masked genome and drop any pair where either
     promoter exceeds a repeat-fraction threshold. Run at several thresholds to
     show the signal is not just transposons.

  3. GC CONTROL. Default is pairwise-drop: discard a pair if the two promoters
     differ in GC by more than the tolerance (0.10, the lab standard). This
     bounds per-comparison GC difference while preserving orthology. We also
     report the raw pool GC gap so you can answer Joel's question ("is GC
     matching even necessary, or are the promoters already close?").

  4. PAIR-LEVEL / CHROMOSOME HOLDOUT. Each ortholog pair is assigned to
     train/val/test by its HUMAN gene's chromosome (EpiModX convention:
     chr8+chr9 -> test, chr10 -> val, rest -> train). Both species' promoters
     of a pair land in the SAME split, so the model can never see a gene's
     human promoter in train and its mouse ortholog in test (ortholog leakage).

Outputs (into --out-dir):
  human_promoters.fasta, <target>_promoters.fasta   # for a FASTA-based harness
  train.csv, dev.csv, test.csv                       # DNABERT-2 native: sequence,label
  manifest.csv                                        # full provenance for auditing
  summary.json                                        # per-stage counts + GC gaps
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass, asdict

# --------------------------------------------------------------------------- #
# Pure, dependency-free helpers (unit-tested in tests/test_core.py)
# --------------------------------------------------------------------------- #

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

TEST_CHROMS = {"chr8", "chr9"}
VAL_CHROMS = {"chr10"}


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def gc_content(seq: str) -> float:
    """GC fraction among A/C/G/T (N ignored in the denominator)."""
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    n = gc + at
    return gc / n if n else 0.0


def repeat_fraction(seq: str) -> float:
    """Fraction of soft-masked (lowercase) bases = RepeatMasker calls."""
    if not seq:
        return 1.0
    return sum(1 for c in seq if c.islower()) / len(seq)


def n_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    return seq.upper().count("N") / len(seq)


def promoter_window(start_1b: int, end_1b: int, strand: str,
                    upstream: int, downstream: int) -> tuple[int, int]:
    """Return a 0-based, half-open [gstart, gend) genomic window around the TSS.

    GTF coords are 1-based inclusive. TSS = gene start on '+', gene end on '-'.
    The returned window, once fetched and (for '-') reverse-complemented, yields
    `upstream` bases 5' of the TSS followed by `downstream` bases into the gene.
    """
    if strand == "+":
        tss0 = start_1b - 1
        gstart = tss0 - upstream
        gend = tss0 + downstream
    elif strand == "-":
        tss0 = end_1b - 1
        gstart = tss0 - downstream + 1
        gend = tss0 + upstream + 1
    else:
        raise ValueError(f"bad strand {strand!r}")
    return max(gstart, 0), gend


def fold_for(human_chrom: str) -> str:
    if human_chrom in TEST_CHROMS:
        return "test"
    if human_chrom in VAL_CHROMS:
        return "val"
    return "train"


def pairwise_gc_keep(pairs: list["Pair"], tol: float) -> list["Pair"]:
    """Keep pairs whose two promoters differ in GC by <= tol."""
    return [p for p in pairs if abs(p.human.gc - p.target.gc) <= tol]


def pool_gc_gap(pairs: list["Pair"]) -> float:
    """Mean(target GC) - mean(human GC) across pairs (signed)."""
    if not pairs:
        return 0.0
    h = sum(p.human.gc for p in pairs) / len(pairs)
    t = sum(p.target.gc for p in pairs) / len(pairs)
    return t - h


# --------------------------------------------------------------------------- #
# Records
# --------------------------------------------------------------------------- #

@dataclass
class Promoter:
    gene_id: str
    chrom: str
    start: int          # 0-based half-open window actually used
    end: int
    strand: str
    seq: str            # oriented 5'->3', case preserved (for repeat calc)
    gc: float
    rep: float

    def model_seq(self) -> str:
        return self.seq.upper()


@dataclass
class Pair:
    human: Promoter
    target: Promoter
    fold: str


# --------------------------------------------------------------------------- #
# IO: GTF parsing + sequence extraction (needs pyfaidx)
# --------------------------------------------------------------------------- #

def parse_gene_tss(gtf_path: str) -> dict[str, tuple[str, int, int, str]]:
    """gene_id (unversioned) -> (chrom, start_1b, end_1b, strand) for gene rows."""
    genes: dict[str, tuple[str, int, int, str]] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            gid = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("gene_id"):
                    gid = field.split('"')[1].split(".")[0].replace("_PAR_Y", "")
                    break
            if gid:
                genes[gid] = (chrom, start, end, strand)
    return genes


def build_promoters(gene_ids, genes, fasta_path, upstream, downstream, max_n):
    """Extract promoter sequences for the requested gene_ids. Returns dict gid->Promoter.

    Windows clipped by a contig end (wrong length) or with too many Ns are skipped.
    """
    import pyfaidx

    fa = pyfaidx.Fasta(fasta_path, sequence_always_upper=False, as_raw=True)
    want = upstream + downstream
    out: dict[str, Promoter] = {}
    for gid in gene_ids:
        meta = genes.get(gid)
        if meta is None:
            continue
        chrom, s1, e1, strand = meta
        if chrom not in fa:
            continue
        gstart, gend = promoter_window(s1, e1, strand, upstream, downstream)
        raw = str(fa[chrom][gstart:gend])
        if len(raw) != want:                      # clipped at contig boundary
            continue
        seq = revcomp(raw) if strand == "-" else raw
        if n_fraction(seq) > max_n:
            continue
        out[gid] = Promoter(gid, chrom, gstart, gend, strand, seq,
                            gc_content(seq), repeat_fraction(seq))
    return out


# --------------------------------------------------------------------------- #
# Assemble
# --------------------------------------------------------------------------- #

def load_orthologs(path: str) -> list[tuple[str, str]]:
    pairs = []
    with open(path) as fh:
        header = fh.readline()  # noqa: F841
        for line in fh:
            h, t = line.rstrip("\n").split("\t")[:2]
            pairs.append((h, t))
    return pairs


def assemble(ortholog_tsv, target_key, out_dir, *, upstream, downstream,
             n_pairs, repeat_max, gc_tol, do_repeat_filter, do_gc_match,
             max_n, seed):
    import os
    from species_config import HUMAN, TARGETS

    target = TARGETS[target_key]
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(seed)
    summary: dict = {"target": target_key, "params": {
        "upstream": upstream, "downstream": downstream, "n_pairs": n_pairs,
        "repeat_max": repeat_max, "gc_tol": gc_tol,
        "repeat_filter": do_repeat_filter, "gc_match": do_gc_match, "seed": seed}}

    ortho = load_orthologs(ortholog_tsv)
    summary["orthologs_1to1"] = len(ortho)

    print("Parsing GTFs ...")
    h_genes = parse_gene_tss(HUMAN.gtf)
    t_genes = parse_gene_tss(target.gtf)

    # keep only pairs where both gene IDs are annotated (with coords) in each GTF
    ortho = [(h, t) for (h, t) in ortho if h in h_genes and t in t_genes]
    summary["orthologs_annotated"] = len(ortho)

    # optional deterministic subsample BEFORE extraction (extraction is the slow part)
    if n_pairs and n_pairs < len(ortho):
        ortho = rng.sample(ortho, n_pairs)
    summary["orthologs_sampled"] = len(ortho)

    print(f"Extracting {len(ortho)} human + {len(ortho)} {target_key} promoters ...")
    h_prom = build_promoters([h for h, _ in ortho], h_genes, HUMAN.genome_fasta,
                             upstream, downstream, max_n)
    t_prom = build_promoters([t for _, t in ortho], t_genes, target.genome_fasta,
                             upstream, downstream, max_n)

    pairs: list[Pair] = []
    for h, t in ortho:
        if h in h_prom and t in t_prom:
            hp, tp = h_prom[h], t_prom[t]
            pairs.append(Pair(hp, tp, fold_for(hp.chrom)))
    summary["pairs_extracted"] = len(pairs)
    summary["raw_gc_gap_target_minus_human"] = round(pool_gc_gap(pairs), 4)

    if do_repeat_filter:
        pairs = [p for p in pairs
                 if p.human.rep <= repeat_max and p.target.rep <= repeat_max]
        summary["pairs_after_repeat_filter"] = len(pairs)

    if do_gc_match:
        pairs = pairwise_gc_keep(pairs, gc_tol)
        summary["pairs_after_gc_match"] = len(pairs)
    summary["gc_gap_final_target_minus_human"] = round(pool_gc_gap(pairs), 4)

    counts = {"train": 0, "val": 0, "test": 0}
    for p in pairs:
        counts[p.fold] += 1
    summary["split_pairs"] = counts
    summary["final_examples_per_class"] = len(pairs)

    _emit(pairs, target_key, out_dir)
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    return summary


def _emit(pairs, target_key, out_dir):
    import csv
    import os

    # FASTA (for a positives/negatives harness) -- explicit distinct names to
    # avoid the classify.py filename-clobbering footgun.
    def _hdr(p, sp):
        return (f">{sp}|{p.gene_id}|{p.chrom}:{p.start}-{p.end}|{p.strand}"
                f"|gc={p.gc:.3f}|rep={p.rep:.3f}")

    with open(os.path.join(out_dir, "human_promoters.fasta"), "w") as fh:
        for p in pairs:
            fh.write(_hdr(p.human, "human") + "\n" + p.human.model_seq() + "\n")
    with open(os.path.join(out_dir, f"{target_key}_promoters.fasta"), "w") as fh:
        for p in pairs:
            fh.write(_hdr(p.target, target_key) + "\n" + p.target.model_seq() + "\n")

    # DNABERT-2 native CSVs (sequence,label). label 0=human, 1=target.
    split_rows = {"train": [], "val": [], "test": []}
    manifest = []
    for p in pairs:
        split_rows[p.fold].append((p.human.model_seq(), 0))
        split_rows[p.fold].append((p.target.model_seq(), 1))
        manifest.append([
            p.fold, p.human.gene_id, p.target.gene_id,
            p.human.chrom, p.human.strand, f"{p.human.gc:.4f}", f"{p.human.rep:.4f}",
            p.target.chrom, p.target.strand, f"{p.target.gc:.4f}", f"{p.target.rep:.4f}",
        ])

    for fold, fname in (("train", "train.csv"), ("val", "dev.csv"), ("test", "test.csv")):
        rows = split_rows[fold]
        with open(os.path.join(out_dir, fname), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sequence", "label"])
            w.writerows(rows)

    with open(os.path.join(out_dir, "manifest.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fold", "human_gene", "target_gene",
                    "h_chrom", "h_strand", "h_gc", "h_rep",
                    "t_chrom", "t_strand", "t_gc", "t_rep"])
        w.writerows(manifest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthologs", required=True, help="TSV from get_orthologs.py")
    ap.add_argument("--target", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--upstream", type=int, default=1000)
    ap.add_argument("--downstream", type=int, default=0)
    ap.add_argument("--n-pairs", type=int, default=2000,
                    help="deterministic subsample of ortholog pairs; 0 = use all")
    ap.add_argument("--repeat-max", type=float, default=0.20,
                    help="drop a pair if either promoter exceeds this repeat fraction")
    ap.add_argument("--gc-tol", type=float, default=0.10,
                    help="drop a pair if |gc_human - gc_target| exceeds this")
    ap.add_argument("--no-repeat-filter", action="store_true")
    ap.add_argument("--no-gc-match", action="store_true")
    ap.add_argument("--max-n", type=float, default=0.0,
                    help="drop a promoter whose N fraction exceeds this")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    assemble(
        args.orthologs, args.target, args.out_dir,
        upstream=args.upstream, downstream=args.downstream,
        n_pairs=(args.n_pairs or None),
        repeat_max=args.repeat_max, gc_tol=args.gc_tol,
        do_repeat_filter=not args.no_repeat_filter,
        do_gc_match=not args.no_gc_match,
        max_n=args.max_n, seed=args.seed,
    )


if __name__ == "__main__":
    main()
