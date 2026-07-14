import re, glob
import numpy as np
import matplotlib.pyplot as plt

# Collect accuracy + per-class AUC from all threeway logs
files = sorted(glob.glob("threeway_seed_*.txt"))
accs, auc_f, auc_p, auc_r = [], [], [], []

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

print(f"Seeds found: {len(accs)}")
print(f"Accuracy: {np.mean(accs):.3f} ± {np.std(accs):.3f}")
print(f"Functional vs rest AUC: {np.mean(auc_f):.3f} ± {np.std(auc_f):.3f}")
print(f"Pseudogene vs rest AUC:  {np.mean(auc_p):.3f} ± {np.std(auc_p):.3f}")
print(f"Random vs rest AUC:      {np.mean(auc_r):.3f} ± {np.std(auc_r):.3f}")

# Bar plot of per-class AUC with error bars
labels = ["Functional\nvs rest", "Pseudogene\nvs rest", "Random\nvs rest"]
means = [np.mean(auc_f), np.mean(auc_p), np.mean(auc_r)]
stds  = [np.std(auc_f), np.std(auc_p), np.std(auc_r)]

plt.figure(figsize=(7,5))
bars = plt.bar(labels, means, yerr=stds, capsize=5,
               color=["#4C72B0", "#C44E52", "#8C8C8C"])
plt.axhline(0.5, color="gray", linestyle="--", alpha=0.6, label="chance (0.5)")
plt.ylabel("One-vs-rest AUC")
plt.title(f"3-way classification ({len(accs)} seeds)\nFunctional vs Pseudogene vs Random promoters")
plt.ylim(0.5, 1.0)
plt.legend()
for b, m in zip(bars, means):
    plt.text(b.get_x()+b.get_width()/2, m+0.01, f"{m:.3f}", ha="center", fontsize=10)
plt.savefig("plot_3way_auc.png", dpi=150, bbox_inches="tight")
print("Saved plot_3way_auc.png")