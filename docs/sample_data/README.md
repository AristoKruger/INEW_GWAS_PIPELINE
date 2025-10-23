# Sample Data Pack

Use these tiny files to practise the pipeline without touching your real project data. Everything lives inside `docs/sample_data/` so you can inspect it safely.

## What is included?
- `raw/raw_genotypes_sample.csv` - three example lines (Sample_1-3) with four SNPs.
- `metadata/35k_probe_set_sample.xlsx` - chromosome and position info for the same SNPs.
- `metadata/35k_probe_set_sample_tab.map` - ready-made PLINK map if you want more PLINK practice.
- `raw/phenotype_sample.txt` - a miniature phenotype table that matches the three example lines.

## Quick practice run
1. Open a terminal in the project folder (`INEW_GWAS_MUDAU_PIPELINE`).
2. Dry run the sample config: `python -m gwas_pipeline --config config/pipeline.sample.yaml --dry-run`
3. Run the pipeline on the sample set (writes to `outputs/sample_run/`): `python -m gwas_pipeline --config config/pipeline.sample.yaml`
4. Open `outputs/run_history/latest.json` to confirm the last finished step, then read the matching `run_*.log` in the same folder.

## Clean up
- Sample outputs live in `outputs/sample_run/`. Delete it whenever you want a fresh practice run.
- The original sample inputs stay in `docs/sample_data/`. Leave them there so the example config keeps working.
