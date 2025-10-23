# step07_extract_labels.py
#
# Extracts genotype labels from a labeled dissimilarity matrix
# Usage: python step07_extract_labels.py <input.csv> <output_labels.txt>
#

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])

# Extract column headers, skip first column (row names)
labels = df.columns[1:]

# Save to file, one label per line
df_out = pd.DataFrame({"Genotype": labels})
df_out.to_csv(sys.argv[2], index=False, header=False)
