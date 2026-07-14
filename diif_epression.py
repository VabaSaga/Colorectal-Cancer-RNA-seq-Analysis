import sqlite3
import pandas as pd
import re
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

folder = r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW"
conn = sqlite3.connect(folder + r"\crc_dex.db")

# --- Step 1: Pull expression data joined with sample metadata ---
query = """
SELECT expression.gene_id, expression.fpkm, samples.patient_id, samples.condition
FROM expression
JOIN samples ON expression.sample_id = samples.sample_id
"""
data = pd.read_sql(query, conn)
conn.close()

# --- Step 2: Remove the corrupted Excel-date gene names ---
date_pattern = re.compile(r'^\d{1,2}-[A-Za-z]{3}$')
data = data[~data["gene_id"].str.match(date_pattern)]

# --- Step 3: Reshape so each row = one gene, one patient, tumor+normal side by side ---
pivoted = data.pivot_table(
    index=["gene_id", "patient_id"],
    columns="condition",
    values="fpkm"
).reset_index()

# Drop any gene/patient rows missing either tumor or normal value
pivoted = pivoted.dropna(subset=["tumor", "normal"])

# --- Step 4: Run a paired Wilcoxon test per gene ---
results = []

for gene_id, group in pivoted.groupby("gene_id"):
    if len(group) < 5:  
        # skip genes with too few paired samples to test meaningfully
        continue
    
    tumor_vals = group["tumor"].values
    normal_vals = group["normal"].values
    
    # Skip genes with no variation at all (test can't run)
    if (tumor_vals == normal_vals).all():
        continue
    
    try:
        stat, p_value = wilcoxon(tumor_vals, normal_vals)
    except ValueError:
        continue
    
    mean_tumor = tumor_vals.mean()
    mean_normal = normal_vals.mean()
    
    results.append({
        "gene_id": gene_id,
        "mean_tumor_fpkm": mean_tumor,
        "mean_normal_fpkm": mean_normal,
        "p_value": p_value
    })

results_df = pd.DataFrame(results)
print(f"Tested {len(results_df)} genes")

# --- Step 5: FDR correction across all tested genes ---
reject, pvals_corrected, _, _ = multipletests(results_df["p_value"], method="fdr_bh")
results_df["padj"] = pvals_corrected
results_df["significant"] = reject

# --- Step 6: Sort by significance and save ---
results_df = results_df.sort_values("padj")

output_path = folder + r"\differential_expression_results.csv"
results_df.to_csv(output_path, index=False)

print(f"Significant genes (padj < 0.05): {results_df['significant'].sum()}")
print(results_df.head(10))
print(f"Saved results to: {output_path}")