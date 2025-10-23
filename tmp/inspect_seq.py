import pandas as pd
from pathlib import Path
path = Path('data/metadata/35k_probe_set_IWGSCv1.xlsx')
df = pd.read_excel(path, usecols=['35K SNPId','Sequence'], nrows=3)
print(df)
for item in df.itertuples(index=False):
    print(item)
