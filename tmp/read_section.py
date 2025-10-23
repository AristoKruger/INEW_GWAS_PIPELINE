from pathlib import Path
text = Path('README.md').read_text(encoding='utf-8')
start = text.index('## Stage 5')
end = text.index('## Stage 6')
print(text[start:end])
