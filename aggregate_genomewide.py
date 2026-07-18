import re
import glob
import statistics

files = sorted(glob.glob("threeway_genomewide_seed_*.txt"))
accs, auc_f, auc_p, auc_r = [], [], [], []
seeds_done = []

for fp in files:
    t = open(fp).read()
    a = re.search(r"Overall accuracy: ([0-9.]+)", t)
    f = re.search(r"AUC functional vs rest: ([0-9.]+)", t)
    p = re.search(r"AUC pseudogene vs rest: ([0-9.]+)", t)
    r = re.search(r"AUC random vs rest: ([0-9.]+)", t)
    if a and f and p and r:
        accs.append(float(a.group(1)))
        auc_f.append(float(f.group(1)))
        auc_p.append(float(p.group(1)))
        auc_r.append(float(r.group(1)))
        seeds_done.append(fp)

print(f"Completed seeds: {len(accs)} ({', '.join(seeds_done)})")
if not accs:
    print("No completed runs yet.")
    raise SystemExit


def summ(name, vals):
    m = statistics.mean(vals)
    s = statistics.stdev(vals) if len(vals) > 1 else 0.0
    print(f"{name:<28} {m:.3f} +/- {s:.3f}  (n={len(vals)})")


summ("Overall accuracy", accs)
summ("AUC functional vs rest", auc_f)
summ("AUC pseudogene vs rest", auc_p)
summ("AUC random vs rest", auc_r)
