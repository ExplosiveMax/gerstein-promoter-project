"""Aggregate the fine-scale chr1 offset scan and plot AUC vs distance from TSS.

Reads chr1_fine_offset_{offset}_seed_{seed}.txt (each has a
  FINAL_RESULT seed=N auc=X auprc=Y
line), aggregates across seeds 0/1/42, and plots the mean curve with a
per-seed spread band and a dashed chance line at 0.5.

Outputs:
  meeting_prep/offset_auc_vs_distance.csv
  meeting_prep/offset_auc_vs_distance.png   (AUROC + AUPRC panels)
"""
import glob
import os
import re
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN_DIR = os.path.join(os.path.dirname(__file__), "task1_offset_run")
SEEDS = [0, 1, 42]
PAT = re.compile(r"FINAL_RESULT seed=(\d+) auc=([\d.]+) auprc=([\d.]+)")

def parse():
    """offset -> {seed: (auc, auprc)}"""
    data = {}
    for path in glob.glob(os.path.join(RUN_DIR, "chr1_fine_offset_*_seed_*.txt")):
        m = re.search(r"offset_(\d+)_seed_(\d+)\.txt", os.path.basename(path))
        if not m:
            continue
        offset = int(m.group(1))
        txt = open(path).read()
        fm = PAT.search(txt)
        if not fm:
            continue
        seed = int(fm.group(1))
        data.setdefault(offset, {})[seed] = (float(fm.group(2)), float(fm.group(3)))
    return data

def main():
    data = parse()
    offsets = sorted(data.keys())
    if not offsets:
        print("No results parsed yet.")
        return

    rows = []
    for off in offsets:
        seed_vals = data[off]
        aucs = [seed_vals[s][0] for s in SEEDS if s in seed_vals]
        aps = [seed_vals[s][1] for s in SEEDS if s in seed_vals]
        rows.append({
            "offset_bp": off,
            "n_seeds": len(aucs),
            "auc_mean": np.mean(aucs), "auc_min": np.min(aucs),
            "auc_max": np.max(aucs), "auc_std": np.std(aucs),
            "auprc_mean": np.mean(aps), "auprc_min": np.min(aps),
            "auprc_max": np.max(aps), "auprc_std": np.std(aps),
            **{f"auc_seed{s}": seed_vals.get(s, (np.nan, np.nan))[0] for s in SEEDS},
            **{f"auprc_seed{s}": seed_vals.get(s, (np.nan, np.nan))[1] for s in SEEDS},
        })

    # ---- CSV ----
    csv_path = os.path.join(os.path.dirname(__file__), "offset_auc_vs_distance.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {csv_path}  ({len(rows)} offsets)")

    x = np.array([r["offset_bp"] for r in rows])
    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
    for ax, metric, color, title in [
        (axes[0], "auc", "#2b6cb0", "AUROC"),
        (axes[1], "auprc", "#c53030", "AUPRC"),
    ]:
        mean = np.array([r[f"{metric}_mean"] for r in rows])
        lo = np.array([r[f"{metric}_min"] for r in rows])
        hi = np.array([r[f"{metric}_max"] for r in rows])
        ax.fill_between(x, lo, hi, color=color, alpha=0.18,
                        label="seed spread (min–max)")
        ax.plot(x, mean, "-o", color=color, lw=2, ms=5,
                label=f"mean {title} (seeds 0/1/42)")
        ax.axhline(0.5, ls="--", color="#4a5568", lw=1.2, label="chance (0.5)")
        ax.set_ylabel(title)
        ax.set_ylim(0.45, max(0.8, hi.max() + 0.03))
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=9)
        ax.set_title(f"{title} vs distance from TSS — chr1, 500bp window, "
                     "promoter vs GC-matched random")
    axes[1].set_xlabel("distance from TSS (bp)")
    fig.tight_layout()
    png = os.path.join(os.path.dirname(__file__), "offset_auc_vs_distance.png")
    fig.savefig(png, dpi=300)
    print(f"wrote {png}")

if __name__ == "__main__":
    main()
