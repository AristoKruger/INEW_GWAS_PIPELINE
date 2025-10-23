from pathlib import Path
text = Path('README.md').read_text(encoding='utf-8')
start = text.index('--dir ')
print(repr(text[start:start+40]))
