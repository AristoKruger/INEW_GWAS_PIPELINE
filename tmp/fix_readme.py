from pathlib import Path
path = Path('README.md')
text = path.read_text(encoding='utf-8')
text = text.replace('outputs\\structure\nuns', 'outputs/structure/runs')
path.write_text(text, encoding='utf-8')
