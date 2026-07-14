import gzip
import glob
import os
import pandas as pd

# Folder where all your GSM files live
folder = r"C:\Users\HP\OneDrive\Desktop\Colorectal Cancer\GSE50760_RAW"

# Find all the .txt.gz files in that folder
files = glob.glob(os.path.join(folder, "*.txt.gz"))
print(f"Found {len(files)} files")

# We'll build a list of small tables, one per sample, then combine them
all_samples = []

for file_path in files:
    # Read the file — it's tab-separated with a header row
    df = pd.read_csv(file_path, sep="\t", compression="gzip")
    
    # The gene column is called "genes", set it as the index
    df = df.set_index("genes")
    
    all_samples.append(df)

# Combine all samples side-by-side, matching on gene name
merged = pd.concat(all_samples, axis=1)

print("Merged shape (genes x samples):", merged.shape)
print(merged.head())

# Save the combined table
output_path = os.path.join(folder, "GSE50760_merged_FPKM.csv")
merged.to_csv(output_path)
print(f"Saved merged file to: {output_path}")

# Keep only tumor (.1) and normal (.2) samples, drop metastasis (.3)
keep_cols = [col for col in merged.columns if ".1_" in col or ".2_" in col]
filtered = merged[keep_cols]

print("Filtered shape (genes x samples):", filtered.shape)
print(filtered.columns.tolist())

# Save this filtered version
filtered_output_path = os.path.join(folder, "GSE50760_tumor_normal_FPKM.csv")
filtered.to_csv(filtered_output_path)
print(f"Saved filtered file to: {filtered_output_path}")

# Build sample metadata: extract patient ID and condition from column names
sample_info = []

for col in filtered.columns:
    # col looks like "AMC_2.1_FPKM" -> patient=2, condition=tumor
    name_part = col.replace("_FPKM", "")   # "AMC_2.1"
    patient_id = name_part.split("_")[1].split(".")[0]   # "2"
    tissue_code = name_part.split(".")[1]                # "1" or "2"
    
    condition = "tumor" if tissue_code == "1" else "normal"
    
    sample_info.append({
        "sample_id": col,
        "patient_id": patient_id,
        "condition": condition
    })

samples_df = pd.DataFrame(sample_info)
print(samples_df)

# Save it
samples_output_path = os.path.join(folder, "samples_metadata.csv")
samples_df.to_csv(samples_output_path, index=False)
print(f"Saved samples metadata to: {samples_output_path}")