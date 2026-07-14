import sqlite3
import pandas as pd
import os

folder = r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW"

# --- Load our two source files ---
expression_wide = pd.read_csv(os.path.join(folder, "GSE50760_tumor_normal_FPKM.csv"))
samples_df = pd.read_csv(os.path.join(folder, "samples_metadata.csv"))

# --- Connect to (and create) the SQLite database ---
db_path = os.path.join(folder, "crc_dex.db")
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# --- Create tables ---
cursor.execute("DROP TABLE IF EXISTS samples")
cursor.execute("DROP TABLE IF EXISTS genes")
cursor.execute("DROP TABLE IF EXISTS expression")

cursor.execute("""
CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    patient_id TEXT,
    condition TEXT
)
""")

cursor.execute("""
CREATE TABLE genes (
    gene_id TEXT PRIMARY KEY
)
""")

cursor.execute("""
CREATE TABLE expression (
    sample_id TEXT,
    gene_id TEXT,
    fpkm REAL,
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
)
""")

conn.commit()
print("Tables created.")

# --- Load samples table ---
samples_df.to_sql("samples", conn, if_exists="append", index=False)
print(f"Loaded {len(samples_df)} rows into samples table.")

# --- Load genes table (deduplicated) ---
genes_df = pd.DataFrame({"gene_id": expression_wide["genes"].unique()})
genes_df.to_sql("genes", conn, if_exists="append", index=False)
print(f"Loaded {len(genes_df)} rows into genes table.")

# --- Reshape expression data from wide to long, then load ---
long_expr = expression_wide.melt(id_vars="genes", var_name="sample_id", value_name="fpkm")
long_expr = long_expr.rename(columns={"genes": "gene_id"})

long_expr.to_sql("expression", conn, if_exists="append", index=False)
print(f"Loaded {len(long_expr)} rows into expression table.")

conn.close()
print("Done. Database saved at:", db_path)