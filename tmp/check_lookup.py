import sys
from pathlib import Path
root = Path(__file__).resolve().parent.parent
sys.path.append(str(root/'workflow'))
sys.path.append(str(root))
import pandas as pd
from workflow.gwas_pipeline.config import PipelineConfig
cfg = PipelineConfig.load(root/'config/pipeline.yaml')
annotation_excel = cfg.path('paths','probe_annotation_excel')
print('Excel path', annotation_excel)
lookup = {}
if annotation_excel and annotation_excel.exists():
    df = pd.read_excel(annotation_excel, usecols=['35K SNPId','Sequence'])
    print('rows', len(df))
    for snp_id, sequence in df.itertuples(index=False):
        if not isinstance(sequence, str):
            continue
        start = sequence.find('[')
        end = sequence.find(']', start+1)
        if start == -1 or end == -1 or end <= start+1:
            continue
        alleles = sequence[start+1:end].split('/')
        if len(alleles) != 2:
            continue
        left = alleles[0].strip().upper()
        right = alleles[1].strip().upper()
        lookup[snp_id] = (left,right)
print('lookup size', len(lookup))
missing = 0
with (root/'outputs/step09_csv_to_tped/filtered_genotypes_strict.tped').open() as fh:
    for i,line in enumerate(fh):
        fields = line.strip().split('\t')
        snp = fields[1]
        if snp not in lookup:
            missing += 1
        if i == 0:
            print('first snp', snp, 'alleles', lookup.get(snp))
print('missing count', missing)
