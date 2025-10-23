from pathlib import Path

text = Path("README.md").read_text(encoding="utf-8")
old = """## Stage 4 ? PCoA and STRUCTURE input

- Run PCoA (R): `python -m gwas_pipeline --steps step13_pcoa`
- Convert to STRUCTURE format: `python -m gwas_pipeline --steps step14_convert_to_structure`

## Stage 5 ? TASSEL GWAS post-processing

- Normalise PED: `python -m gwas_pipeline --steps step15_fix_pruned_ped`
- Split traits: `python -m gwas_pipeline --steps step16_trait_split`
- Adjust p-values and summarise: `python -m gwas_pipeline --steps step17_trait_adjust`

After step17 a short report is written to `outputs/summary/report.md`.

"""
print(old in text)
