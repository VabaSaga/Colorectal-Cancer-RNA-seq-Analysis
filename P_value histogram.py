import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Load data ---
res_df = pd.read_csv(r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW\differential_expression_results.csv", index_col=0)

# --- 2. Plot histogram of raw p-values ---
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(res_df["p_value"].dropna(), bins=50, color="steelblue", edgecolor="white")
ax.set_xlabel("p-value")
ax.set_ylabel("Frequency")
ax.set_title("P-value Distribution")
plt.tight_layout()
plt.savefig("pvalue_histogram.png", dpi=300, bbox_inches="tight")
plt.show()