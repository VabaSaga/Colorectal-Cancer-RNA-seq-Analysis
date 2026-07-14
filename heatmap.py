import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# --- 1. Load data ---
res_df = pd.read_csv(r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW\differential_expression_results.csv", index_col=0)

# --- 2. Compute log2 fold change (tumor vs normal) ---
pseudo = 0.01  # avoids log(0)
res_df["log2FC"] = np.log2(
    (res_df["mean_tumor_fpkm"] + pseudo) / (res_df["mean_normal_fpkm"] + pseudo)
)

# --- 3. Select top genes among significant ones, ranked by |log2FC| ---
sig = res_df[res_df["significant"] == True].dropna(subset=["padj"])
top_genes = sig.reindex(sig["log2FC"].abs().sort_values(ascending=False).index).head(20)

# --- 4. Sort by log2FC for display (descending) ---
top_genes = top_genes.sort_values("log2FC", ascending=False)

# --- 5. Single-column heatmap colored by log2FC ---
lfc_col = top_genes[["log2FC"]]

fig, ax = plt.subplots(figsize=(3, 8))
sns.heatmap(
    lfc_col,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "log2FC (Tumor vs Normal)"},
    yticklabels=True,
    ax=ax
)
ax.set_title("Top 20 DE Genes\nby log2FC", fontsize=11)
ax.set_ylabel("")
ax.set_xticklabels(["log2FC"])
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig("de_heatmap_log2fc.png", dpi=300, bbox_inches="tight")
plt.show()