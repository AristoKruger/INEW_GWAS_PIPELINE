import pandas as pd
from pathlib import Path
path = Path('data/metadata/35k_probe_set_IWGSCv1.xlsx')
row = pd.read_excel(path, usecols=['35K SNPId','Sequence'], nrows=1, skiprows=lambda x: x not in [0])
print(row)
