import pandas as pd
from pathlib import Path
path = Path('data/metadata/35k_probe_set_IWGSCv1.xlsx')
df = pd.read_excel(path, usecols=['35K SNPId','Sequence'])
row = df.loc[df['35K SNPId']=='AX-94977015'].head(1)
print(row.to_string(index=False))
