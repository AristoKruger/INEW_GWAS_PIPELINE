import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python step04_transpose_matrix.py <input_csv> <output_csv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the genotype file with headers
df = pd.read_csv(input_file)

# Set 'Probe ID' as the index so it becomes column headers after transpose
df.set_index('Probe ID', inplace=True)

# Transpose the matrix
df_T = df.transpose()

# Reset index to make taxa a column
df_T.reset_index(inplace=True)
df_T.rename(columns={"index": "Taxa"}, inplace=True)

# Save transposed file
df_T.to_csv(output_file, index=False)

print(f"✅ Transposed genotype file saved to {output_file}")

