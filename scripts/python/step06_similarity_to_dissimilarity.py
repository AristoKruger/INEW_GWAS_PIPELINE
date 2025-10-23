# 09_similarity_to_dissimilarity.py
#
# Converts a pairwise similarity matrix to a pairwise dissimilarity matrix (1 - similarity)
# Usage: python 09_similarity_to_dissimilarity.py <similarity_matrix.csv> <output_dissimilarity.csv>
#

import sys
import pandas as pd

# Read in similarity matrix
similarities = pd.read_csv(sys.argv[1])

# Copy the matrix for modification
dissimilarities = similarities.copy()

# Get list of genotype columns (excluding the first column which is row labels)
genotypes = dissimilarities.columns[1:]

# Subtract from 1 to get dissimilarity values
dissimilarities[genotypes] = 1 - dissimilarities[genotypes]

# Save to output
dissimilarities.to_csv(sys.argv[2], index=False)
