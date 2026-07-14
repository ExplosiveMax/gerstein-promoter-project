import os
import re
import statistics
import matplotlib.pyplot as plt

def get_auc(filepath):
    """Pull the AUC out of a log file via the FINAL_RESULT line."""
    with open(filepath) as f:
        text = f.read()
    m = re.search(r"FINAL_RESULT.*auc=([0-9.]+)", text)
    return float(m.group(1)) if m else None

def collect(folder, prefix, values, seeds=(0, 1, 42, 123, 999)):
    """For each value (window size or offset), average AUC across seeds."""
    means, stds = [], []
    for v in values:
        aucs = []
        for s in seeds:
            path = os.path.join(folder, f"{prefix}_{v}_seed_{s}.txt")
            if os.path.exists(path):
                auc = get_auc(path)
                if auc is not None:
                    aucs.append(auc)
        means.append(statistics.mean(aucs) if aucs else None)
        stds.append(statistics.stdev(aucs) if len(aucs) > 1 else 0)
        print(f"{prefix} {v}: mean={means[-1]:.4f} (n={len(aucs)})")
    return means, stds

# ---- Window size plot ----
window_sizes = [256, 500, 1000, 2000]
print("=== Window size ===")
w_means, w_stds = collect("window_output", "window", window_sizes)

plt.figure(figsize=(7, 5))
plt.errorbar(window_sizes, w_means, yerr=w_stds, marker="o", capsize=4)
plt.xlabel("Region size (bp)")
plt.ylabel("Mean AUC (5 seeds)")
plt.title("Promoter vs Random: AUC vs Region Size")
plt.grid(True, alpha=0.3)
plt.ylim(0.5, 1.0)
plt.savefig("plot_window_size.png", dpi=150, bbox_inches="tight")
print("Saved plot_window_size.png")

# ---- Offset plot ----
offsets = [-1000, -500, 0, 500, 1000]
print("\n=== Offset ===")
o_means, o_stds = collect("offset_output", "offset", offsets)

plt.figure(figsize=(7, 5))
plt.errorbar(offsets, o_means, yerr=o_stds, marker="o", capsize=4, color="darkorange")
plt.axvline(0, color="gray", linestyle="--", alpha=0.5, label="TSS")
plt.xlabel("Window offset relative to TSS (bp)")
plt.ylabel("Mean AUC (5 seeds)")
plt.title("Promoter vs Random: AUC vs Window Position")
plt.grid(True, alpha=0.3)
plt.ylim(0.5, 1.05)
plt.legend()
plt.savefig("plot_offset.png", dpi=150, bbox_inches="tight")
print("Saved plot_offset.png")