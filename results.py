import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load your actual data
# Update 'differential_expression_results.csv' if the path to your file is different
df = pd.read_csv("differential_expression_results.csv")

# 2. Calculate log2 Fold Change (adding a small epsilon to handle divisions by zero)
epsilon = 1e-5
df['log2FC'] = np.log2((df['mean_tumor_fpkm'] + epsilon) / (df['mean_normal_fpkm'] + epsilon))

# Calculate -log10 of adjusted p-value (padj is more robust for publication than p_value)
df['-log10_padj'] = -np.log10(df['padj'])

# 3. Establish Thresholds for Highlighting
FC_threshold = 1.0       # Corresponds to a 2-fold change
padj_threshold = 0.05    # Alpha level for statistical significance

# Categorize genes for coloring
df['Expression'] = 'Not Significant'
df.loc[(df['log2FC'] > FC_threshold) & (df['padj'] < padj_threshold), 'Expression'] = 'Upregulated in Tumor'
df.loc[(df['log2FC'] < -FC_threshold) & (df['padj'] < padj_threshold), 'Expression'] = 'Downregulated in Tumor'

# 4. Create Plot Architecture
plt.figure(figsize=(9, 7), dpi=300) # 300 DPI is standard publication quality
sns.set_style("ticks")

# Custom publication-grade color palette
colors = {
    'Not Significant': '#bcbcbc',       # Soft grey
    'Upregulated in Tumor': '#e41a1c',   # Deep Red
    'Downregulated in Tumor': '#377eb8'  # Deep Blue
}

# Scatter plot
sns.scatterplot(
    data=df,
    x='log2FC',
    y='-log10_padj',
    hue='Expression',
    palette=colors,
    alpha=0.75,
    edgecolor=None,
    s=20
)

# Add threshold guideline cutoffs
plt.axvline(x=FC_threshold, color='dimgray', linestyle='--', linewidth=1, alpha=0.8)
plt.axvline(x=-FC_threshold, color='dimgray', linestyle='--', linewidth=1, alpha=0.8)
plt.axhline(y=-np.log10(padj_threshold), color='dimgray', linestyle='--', linewidth=1, alpha=0.8)

# 5. Automatically label top significant genes (Top 5 up, Top 5 down)
# Filters out infinitely high -log10 values if padj was 0
finite_df = df[np.isfinite(df['-log10_padj'])]
top_up = finite_df[finite_df['Expression'] == 'Upregulated in Tumor'].nsmallest(5, 'padj')
top_down = finite_df[finite_df['Expression'] == 'Downregulated in Tumor'].nsmallest(5, 'padj')

for _, row in pd.concat([top_up, top_down]).iterrows():
    plt.text(
        x=row['log2FC'],
        y=row['-log10_padj'] + 0.2, # slight vertical offset
        s=row['gene_id'],
        fontsize=8,
        weight='semibold',
        ha='center',
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1)
    )

# 6. Polishing Labels & Aesthetics
plt.title("Tumor vs. Normal Differential Gene Expression", fontsize=14, pad=15, weight='bold')
plt.xlabel(r"$\log_2 \text{ Fold Change (Tumor / Normal)}$", fontsize=12, labelpad=8)
plt.ylabel(r"$-\log_{10}(\text{Adjusted p-value})$", fontsize=12, labelpad=8)

# Clean up legend and spine layout
plt.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none')
sns.despine()
plt.tight_layout()

# Save for GitHub/Manuscript
plt.savefig("volcano_plot_publication.png", dpi=300, bbox_inches='tight')
plt.show()