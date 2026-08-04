"""Summarize the macaque cross-species run: decomposition + holdout.

Reads:
  runs/macaque/{raw,repeat,repeat_gc}/seed_{0,1,42}.txt   (FINAL_RESULT ... auc= auprc=)
  runs/macaque/repeat_gc/holdout_seed_{0,1,42}.txt        (HOLDOUT_RESULT ... auroc= auprc=)
  runs/macaque/{level}/summary.json                       (pair counts)

Prints a table of mean AUROC + AUPRC across seeds per level, plus the holdout
numbers, and writes it to runs/macaque/SUMMARY.txt.
"""
import glob
import json
import os
import re

RUN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "macaque")
LEVELS = ["raw", "repeat", "repeat_gc"]
SEEDS = [0, 1, 42]

FINAL = re.compile(r"FINAL_RESULT\s+seed=(\d+)\s+auc=([\d.]+)(?:\s+auprc=([\d.]+))?")
HOLD = re.compile(r"HOLDOUT_RESULT\s+dir=\S+\s+seed=(\d+)\s+auroc=([\d.]+)\s+auprc=([\d.]+)")


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def spread(xs):
    return (min(xs), max(xs)) if xs else (float("nan"), float("nan"))


def parse_decomp():
    out = {}
    for level in LEVELS:
        per = {}
        for seed in SEEDS:
            p = os.path.join(RUN, level, f"seed_{seed}.txt")
            if not os.path.exists(p):
                continue
            m = FINAL.search(open(p).read())
            if m:
                auc = float(m.group(2))
                ap = float(m.group(3)) if m.group(3) else None
                per[seed] = (auc, ap)
        out[level] = per
    return out


def parse_holdout():
    per = {}
    for seed in SEEDS:
        p = os.path.join(RUN, "repeat_gc", f"holdout_seed_{seed}.txt")
        if not os.path.exists(p):
            continue
        m = HOLD.search(open(p).read())
        if m:
            per[seed] = (float(m.group(2)), float(m.group(3)))
    return per


def pair_counts():
    counts = {}
    for level in LEVELS:
        p = os.path.join(RUN, level, "summary.json")
        if os.path.exists(p):
            d = json.load(open(p))
            counts[level] = d.get("final_examples_per_class")
    return counts


def fmt_mean_spread(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "   n/a "
    lo, hi = spread(vals)
    return f"{mean(vals):.4f} [{lo:.4f}-{hi:.4f}]"


def main():
    dec = parse_decomp()
    hold = parse_holdout()
    counts = pair_counts()

    lines = []
    lines.append("=" * 78)
    lines.append("MACAQUE cross-species promoter classification — human vs rhesus (Mmul_10)")
    lines.append("n_pairs=1000 requested; seeds 0/1/42; MPS (non-deterministic → spread shown)")
    lines.append("=" * 78)
    lines.append("")
    lines.append("A) Confound decomposition  (random 80/20 split, run_experiment.sh)")
    lines.append(f"   {'level':<10} {'pairs/class':>11}   {'mean AUROC [min-max]':<24}  {'mean AUPRC [min-max]':<24}")
    for level in LEVELS:
        per = dec.get(level, {})
        aucs = [per[s][0] for s in SEEDS if s in per]
        aps = [per[s][1] for s in SEEDS if s in per]
        n = counts.get(level, "?")
        lines.append(f"   {level:<10} {str(n):>11}   {fmt_mean_spread(aucs):<24}  {fmt_mean_spread(aps):<24}")
    lines.append("")
    lines.append("   per-seed AUROC/AUPRC:")
    for level in LEVELS:
        per = dec.get(level, {})
        cells = "  ".join(
            f"s{s}={per[s][0]:.3f}/{(f'{per[s][1]:.3f}' if per[s][1] is not None else 'na')}"
            for s in SEEDS if s in per
        )
        lines.append(f"     {level:<10} {cells}")
    lines.append("")
    lines.append("B) Chromosome-holdout eval on repeat_gc  (classify_holdout.py; chr8/9 test)")
    aucs = [hold[s][0] for s in SEEDS if s in hold]
    aps = [hold[s][1] for s in SEEDS if s in hold]
    lines.append(f"   mean AUROC {fmt_mean_spread(aucs)}   mean AUPRC {fmt_mean_spread(aps)}")
    for s in SEEDS:
        if s in hold:
            lines.append(f"     seed {s:<2}  auroc={hold[s][0]:.4f}  auprc={hold[s][1]:.4f}")
    lines.append("")
    lines.append("Context (prior species, repeat_gc holdout): mouse ~0.94, chimp ~0.50.")
    lines.append("Macaque (~25 My) is expected to fall between them.")
    lines.append("=" * 78)

    text = "\n".join(lines)
    print(text)
    with open(os.path.join(RUN, "SUMMARY.txt"), "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
