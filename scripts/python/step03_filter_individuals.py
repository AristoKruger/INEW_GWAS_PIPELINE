import pandas as pd
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python step03_filter_individuals.py <input_genotype_csv> <output_csv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Derive log file path
log_file = os.path.splitext(output_file)[0] + "_removed_individuals.txt"

# Load genotype file
df = pd.read_csv(input_file)

# Separate probe IDs and genotype matrix
probe_ids = df.iloc[:, 0]
genotypes = df.iloc[:, 1:]

# Define missing code
missing_code = 9

# Compute missing rate per individual
missing_percent = (genotypes == missing_code).sum() / genotypes.shape[0]

# Identify individuals to retain and remove
columns_to_keep = missing_percent[missing_percent < 0.05].index.tolist()
columns_removed = missing_percent[missing_percent >= 0.05].index.tolist()

# Filter genotype matrix
filtered_df = pd.concat([probe_ids, genotypes[columns_to_keep]], axis=1)

# Save filtered genotype data
filtered_df.to_csv(output_file, index=False)

# Write removed individuals to log
with open(log_file, "w") as f:
    f.write("Individuals removed due to >= 5% missing data:\n")
    for name in columns_removed:
        pct = round(missing_percent[name] * 100, 2)
        f.write(f"{name}\t{pct}% missing\n")

# Print summary
print(f"✅ Total individuals before filtering: {genotypes.shape[1]}")
print(f"✅ Total individuals after filtering:  {len(columns_to_keep)}")
print(f"📁 Filtered file saved to: {output_file}")
print(f"🗑️  Log of removed individuals saved to: {log_file}")

