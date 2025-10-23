from pathlib import Path
path = Path('README.md')
text = path.read_text(encoding='utf-8')
text = text.replace('--recodeA `\n      --tab', '--recode `\n      --tab\n      --alleleACGT')
path.write_text(text, encoding='utf-8')
