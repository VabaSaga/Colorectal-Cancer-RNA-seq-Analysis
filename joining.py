import sqlite3
conn = sqlite3.connect(r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW\crc_dex.db")

import pandas as pd

query = """
SELECT expression.gene_id, samples.condition, AVG(expression.fpkm) as avg_fpkm
FROM expression
JOIN samples ON expression.sample_id = samples.sample_id
GROUP BY expression.gene_id, samples.condition
"""

result = pd.read_sql(query, conn)
print(result)
import re

# Get all unique gene names from the genes table
genes_query = "SELECT gene_id FROM genes"
all_genes = pd.read_sql(genes_query, conn)

# Pattern to catch things like "01-Dec", "01-Mar", "01-Sep", "2-Mar", etc.
date_pattern = re.compile(r'^\d{1,2}-[A-Za-z]{3}$')

suspicious = all_genes[all_genes["gene_id"].str.match(date_pattern)]

print(f"Found {len(suspicious)} suspicious gene names out of {len(all_genes)} total")
print(suspicious["gene_id"].unique())

# Remove suspicious date-like gene names from further analysis
clean_genes = all_genes[~all_genes["gene_id"].str.match(date_pattern)]

print(f"Kept {len(clean_genes)} clean gene names, removed {len(all_genes) - len(clean_genes)}")

# Save the list of excluded genes for documentation purposes
suspicious["gene_id"].to_csv(
    r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW\excluded_genes_excel_bug.csv",
    index=False
)


