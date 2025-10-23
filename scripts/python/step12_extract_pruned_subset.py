# step12_extract_pruned_subset.py
#
# Extracts subset of SNPs based on PLINK LD prune output list (.prune.in)
# Usage: python step12_extract_pruned_subset.py <.prune.in> <input.csv> <out.csv>
#

import sys

prune_in_file = sys.argv[1]       # e.g., ld_pruned_genotypes_strict.prune.in
input_csv = sys.argv[2]           # e.g., filtered_genotypes_strict.csv
output_csv = sys.argv[3]          # e.g., ld_pruned_genotypes_strict.csv

# Load selected SNPs into a set for fast lookup
with open(prune_in_file, "r") as f:
    selected_snps = set(line.strip() for line in f)

# Filter rows from the CSV
with open(input_csv, "r") as fin, open(output_csv, "w") as fout:
    for i, line in enumerate(fin):
        if i == 0 or line.split(",")[0] in selected_snps:
            fout.write(line)
