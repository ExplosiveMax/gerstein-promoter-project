"""Report AUROC *and* AUPRC (+ a bootstrap CI) from a predictions file.

Predictions file: CSV with columns `label,score` where score is P(class==1),
i.e. P(target species). Balanced classes here, so AUPRC baseline is ~0.5.

    python eval_metrics.py preds.csv
    python eval_metrics.py preds.csv --n-boot 2000
"""

import argparse
import csv
import random

from sklearn.metrics import roc_auc_score, average_precision_score


def load(path):
    y, s = [], []
    with open(path) as fh:
        r = csv.DictReader(fh)
        for row in r:
            y.append(int(row["label"]))
            s.append(float(row["score"]))
    return y, s


def boot_ci(y, s, metric, n_boot, seed=0):
    rng = random.Random(seed)
    n = len(y)
    vals = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        yb = [y[i] for i in idx]
        if len(set(yb)) < 2:
            continue
        vals.append(metric(yb, [s[i] for i in idx]))
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[int(0.975 * len(vals))]
    return lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("preds")
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()

    y, s = load(args.preds)
    auroc = roc_auc_score(y, s)
    auprc = average_precision_score(y, s)
    rl, rh = boot_ci(y, s, roc_auc_score, args.n_boot)
    pl, ph = boot_ci(y, s, average_precision_score, args.n_boot)
    print(f"n={len(y)}  pos={sum(y)}")
    print(f"AUROC {auroc:.4f}  95% CI [{rl:.4f}, {rh:.4f}]")
    print(f"AUPRC {auprc:.4f}  95% CI [{pl:.4f}, {ph:.4f}]  (baseline ~{sum(y)/len(y):.3f})")


if __name__ == "__main__":
    main()
