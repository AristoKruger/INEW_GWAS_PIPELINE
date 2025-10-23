from pathlib import Path
import pandas as pd
path = Path('outputs/step11_pruned_dataset/ld_pruned_genotypes_strict.csv')
df = pd.read_csv(path)
print(df.columns[:5])
