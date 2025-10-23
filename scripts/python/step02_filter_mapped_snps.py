import pandas as pd
import sys

if len(sys.argv) != 4:
    print("Usage: python step02_filter_mapped_snps.py <input_genotype_csv> <iwgsc_probe_excel> <output_csv>")
    sys.exit(1)

input_file = sys.argv[1]
probe_file = sys.argv[2]
output_file = sys.argv[3]

# Load genotype file
geno_df = pd.read_csv(input_file)

# Load probe annotation file
probe_df = pd.read_excel(probe_file)

# Get list of SNPs with valid chromosome assignments
mapped_snps = probe_df[probe_df['IWGSC_v1_Chromosome'].notna()]['35K SNPId'].unique()

# Filter the genotype dataframe
filtered_geno_df = geno_df[geno_df['Probe ID'].isin(mapped_snps)]

# Output result
filtered_geno_df.to_csv(output_file, index=False)

print(f"✅ SNPs before filtering: {geno_df.shape[0]}")
print(f"✅ SNPs after filtering:  {filtered_geno_df.shape[0]}")
print(f"📁 Filtered genotype file saved to {output_file}")

