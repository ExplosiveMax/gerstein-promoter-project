import re, glob
import numpy as np
import matplotlib.pyplot as plt

# Collect accuracy + per-class AUC from all genome-wide threeway logs
files = sorted(glob.glob("threeway_genomewide_seed_*.txt"))
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
print(f"Accuracy: {np.mean(accs):.3f} +/- {np.std(accs):.3f}")
print(f"Functional vs rest AUC: {np.mean(auc_f):.3f} +/- {np.std(auc_f):.3f}")
print(f"Pseudogene vs rest AUC:  {np.mean(auc_p):.3f} +/- {np.std(auc_p):.3f}")
print(f"Random vs rest AUC:      {np.mean(auc_r):.3f} +/- {np.std(auc_r):.3f}")

# Bar plot: genome-wide vs chr22-only per-class AUC, side by side.
# chr22-only means are the published 5-seed numbers from threeway_seed_*.txt.
chr22_means = [0.867, 0.850, 0.924]   # functional, pseudogene, random (vs rest)
gw_means = [np.mean(auc_f), np.mean(auc_p), np.mean(auc_r)]
gw_stds  = [np.std(auc_f), np.std(auc_p), np.std(auc_r)]

labels = ["Functional\nvs rest", "Pseudogene\nvs rest", "Random\nvs rest"]
x = np.arange(len(labels))
w = 0.38

plt.figure(figsize=(8, 5))
plt.bar(x - w/2, chr22_means, w, label="chr22 only (gene_type match)",
        color="#B0B0B0")
plt.bar(x + w/2, gw_means, w, yerr=gw_stds, capsize=5,
        label="genome-wide (Yale-UCSC 2-way consensus)", color="#4C72B0")
plt.axhline(0.5, color="gray", linestyle="--", alpha=0.6, label="chance (0.5)")
plt.xticks(x, labels)
plt.ylabel("One-vs-rest AUC")
plt.title(f"3-way promoter classification: chr22 vs genome-wide\n"
          f"(genome-wide n={len(accs)} seeds)")
plt.ylim(0.5, 1.0)
plt.legend(fontsize=8)
for xi, m in zip(x + w/2, gw_means):
    plt.text(xi, m + 0.01, f"{m:.3f}", ha="center", fontsize=9)
plt.savefig("plot_3way_genomewide_auc.png", dpi=150, bbox_inches="tight")
print("Saved plot_3way_genomewide_auc.png")
