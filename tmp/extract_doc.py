import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

docx_path = Path('docs/SU-PBL GWAS Pipeline.docx')
with zipfile.ZipFile(docx_path) as z:
    xml = z.read('word/document.xml')
root = ET.fromstring(xml)
ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
texts = []
for para in root.findall('.//w:p', ns):
    runs = [t.text for t in para.findall('.//w:t', ns) if t.text]
    if runs:
        texts.append(''.join(runs))
Path('tmp').mkdir(exist_ok=True)
with open('tmp/_doc_extract.txt', 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(texts))
