import pandas as pd
from pathlib import Path
path = Path('data/metadata/35k_probe_set_IWGSCv1.xlsx')
print(path.exists())
if path.exists():
    df = pd.read_excel(path, nrows=5)
    print(df.columns.tolist())
    print(df.head().to_string())
