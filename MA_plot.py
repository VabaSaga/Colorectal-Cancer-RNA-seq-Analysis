import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Load data ---
res_df = pd.read_csv(r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW\differential_expression_results.csv", index_col=0)

# --- 2. Compute mean expression and log2FC ---
pseudo = 0.01
res_df["mean_expr"] = (res_df["mean_tumor_fpkm"] + res_df["mean_normal_fpkm"]) / 2
res_df["log2FC"] = np.log2((res_df["mean_tumor_fpkm"] + pseudo) / (res_df["mean_normal_fpkm"] + pseudo))

# --- 3. Plot ---
fig, ax = plt.subplots(figsize=(7, 5))

# Non-significant genes (background, gray)
non_sig = res_df[res_df["significant"] == False]
ax.scatter(non_sig["mean_expr"], non_sig["log2FC"],
           s=5, color="gray", alpha=0.4, label="Not significant")

# Significant genes (highlighted, red)
sig = res_df[res_df["significant"] == True]
ax.scatter(sig["mean_expr"], sig["log2FC"],
           s=5, color="red", alpha=0.6, label="Significant (padj < 0.05)")

ax.set_xscale("log")
ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
ax.set_xlabel("Mean Expression (FPKM, log scale)")
ax.set_ylabel("log2 Fold Change (Tumor vs Normal)")
ax.set_title("MA Plot")
ax.legend(markerscale=3, fontsize=9)

plt.tight_layout()
plt.savefig("MA_plot.png", dpi=300, bbox_inches="tight")
plt.show()