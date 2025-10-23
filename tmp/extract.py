from pathlib import Path
text = Path('README.md').read_text(encoding='utf-8')
start = text.index('plink --bfile outputs/step10_pruned_dataset/pruned_panel')
print(text[start:start+300])
