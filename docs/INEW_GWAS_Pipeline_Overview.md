# INEW GWAS Pipeline Overview

_Comprehensive summary of the workflow implementation, configuration, directories, inputs, outputs, and per-step context._

---

## 1. Executive Summary

The INEW GWAS pipeline is a reproducible, hybrid Python/R/PLINK workflow that processes raw SNP microarray genotypes through quality control, LD pruning, TASSEL-based GWAS, optional STRUCTURE analysis, and post-GWAS reporting. The command `python -m gwas_pipeline` orchestrates automated stages defined in `config/pipeline.yaml`, while key manual checkpoints (PLINK helper, STRUCTURE, TASSEL) rely on documented directories and conventions.

---

## 2. Directory Structure

| Path | Purpose |
| --- | --- |
| `config/` | YAML configurations (primary `pipeline.yaml`, sample/demo configs). |
| `data/raw/` | Input genotype CSVs and phenotype tables (read-only after validation). |
| `data/metadata/` | Probe annotation Excel (`35k_probe_set_IWGSCv1.xlsx`) and cleaned map (`35k_probe_set_IWGSCv1_cleaned_tab.map`). |
| `docs/` | User-facing documents (Word specification, presentation flowchart, sample data). |
| `gwas_pipeline/` | Lightweight package wrapper enabling `python -m gwas_pipeline`. |
| `workflow/gwas_pipeline/` | Core Python implementation (`pipeline.py`, `steps.py`, `config.py`, `progress.py`, `utils.py`). |
| `scripts/python/` | Helper scripts (`helper_self_check.py`, validation utilities). |
| `scripts/r/` | R scripts for LD decay (`step10_ld_decay_enhanced.R`), PCoA, TASSEL post-processing. |
| `scripts/cli_examples/` | PowerShell helpers (notably `run_plink_stage.ps1`). |
| `outputs/stepNN_<name>/` | Per-step artefacts (numbered directories). |
| `outputs/run_history/` | Logs and JSON progress snapshots per execution. |
| `outputs/structure/` | STRUCTURE workspace: `raw/`, `runs/<DATE>_K#/`, `structure_harvester/harvester_input/`, `structure_harvester/harvester_output/`, helper scripts under `structure_harvester/structureHarvester_scripts/`, and CLUMPAK `input/`+`outputs/`. |
| `outputs/tassel/` | TASSEL workspace: `diagnostics/` (genotype summaries, kinship matrices, MDS outputs, `mlm_input_table.txt`), `phenotypes/`, `plots/`, `derived/`, `project/`, plus top-level `mlm_stats.txt`, `mlm_effects.txt`, optional `_README.txt`. |

---

## 3. Naming Conventions & Artefacts

- **Step directories**: `outputs/stepNN_<slug>/` align with configuration keys (`step01_filter_variants`, -, `step17_trait_adjust`).
- **Manual outputs**:
  - PLINK helper (Stage 3) writes `pruned_panel.*`, `.prune.in`, and logs under `outputs/step11_plink_pruning/`, then rerun Steps 04-07/13 with pruned inputs to populate `outputs/step13_pcoa/pcoa_coordinates_pruned.csv` / `pcoa_plot_pruned.png`.
  - STRUCTURE stage stores raw runs in `outputs/structure/runs/`, HARVESTER files inside `structure_harvester/{harvester_input,harvester_output}`, helper scripts under `structure_harvester/structureHarvester_scripts/`, and CLUMPAK bundles in `clumpak/input/` (zipped runs) and `clumpak/outputs/` (consensus plots).
  - TASSEL stage exports diagnostics under `outputs/tassel/diagnostics/`, keeps phenotypes in `outputs/tassel/phenotypes/`, plots in `outputs/tassel/plots/`, MLM tables in `outputs/tassel/mlm_stats.txt` / `mlm_effects.txt`, and you may maintain `_README.txt` for run notes.
- **Logs**: `outputs/run_history/run_YYYYMMDD_HHMMSS[_dryrun].log/json`; `latest.json` references the most recent run.
- **Configs**: Each step block in `pipeline.yaml` contains explicit `input`, `output`, numeric parameters, and tool-specific settings; the optional `pipeline.pruned.yaml` adjusts paths for pruned-only iterations.

---

## 4. Inputs & Prerequisites

| Input | Location | Notes |
| --- | --- | --- |
| Genotype matrix | `data/raw/genotypes/raw_genotypes.csv` | SNP rows with `Probe ID` and 0/1/2/9 codes. |
| Probe annotation (Excel) | `data/metadata/35k_probe_set_IWGSCv1.xlsx` | Supplies chromosome mapping for Step 02. |
| Cleaned probe map | `data/metadata/35k_probe_set_IWGSCv1_cleaned_tab.map` | Used to generate PLINK `.map`. |
| Phenotypes for TASSEL | `outputs/tassel/phenotypes/<run_label>.txt` | Manual preparation; TASSEL import. |
| External executables | PLINK, Rscript | Must be accessible via PATH or environment variables (`PLINK_PATH`, `RSCRIPT_PATH`). |

Helper: `python scripts/python/helper_self_check.py` confirms the presence of critical inputs before running the pipeline.

---

## 5. Pipeline Configuration Highlights (`config/pipeline.yaml`)

- `project_root`: relative path resolution for outputs.
- `paths`: canonical locations for raw genotypes, annotations, TASSEL MLM stats, LD summaries, prune lists.
- `steps.stepNN_*`: individual dictionaries specifying input/output paths, thresholds (MAF, missingness), script invocations, etc.
- `tools`: PLINK and LD decay parameter defaults (`plink_ld_prune`, `plink_ld_decay`) plus optional executable overrides.
- Steps mirror the order defined in `workflow/gwas_pipeline/pipeline.py` (`STEP_SEQUENCE`).

Dry-run command `python -m gwas_pipeline --dry-run` validates configuration without file generation.

---

## 6. Scripts & Automation

| Script | Location | Purpose |
| --- | --- | --- |
| `helper_self_check.py` | `scripts/python/` | Validates required inputs and previously generated artefacts. |
| `run_plink_stage.ps1` | `scripts/cli_examples/` | Executes PLINK pruning sequence (indep-pairwise, make-bed, recode). |
| `step10_ld_decay_enhanced.R` | `scripts/r/` | LD-decay helper sourced by Step 10 to summarise distance vs r² and render the decay plot. |
| `step13_pcoa.R` | `scripts/r/` | Performs classical MDS/PCoA on dissimilarity matrix. |
| `step16_trait_split.R` / `step17_trait_adjust.R` | `scripts/r/` | Post-TASSEL processing (per-trait split, multiple-testing adjustments). |
| `workflow/gwas_pipeline/steps.py` | - | Houses all Python step implementations, R/PLINK shellouts, logging, and summary utilities. |

All automated steps use `step_logger` for timing and rely on `PipelineConfig` helpers for path resolution.

---

## 7. Pipeline Stages & Context

### Stage 0 - Workspace Setup
- Create virtual environment (`python -m venv .venv`; activate).
- Install requirements (`pip install -r workflow/requirements.txt`).
- Ensure PLINK and Rscript are accessible; set `PLINK_PATH`/`Path` variables if necessary.
- Install required R packages (`dplyr`, `stringr`, `ggplot2`).
- Run `helper_self_check.py`.

### Stage 1 - QC & Formatting (`step01`-`step05`)
- **Inputs**: raw genotype CSV, annotation Excel/map.
- **Actions**: SNP MAF/missing filtering, map-based retention, individual filtering, matrix transpose, IBS similarity.
- **Outputs**: `filtered_genotypes_strict.csv`, `similarity_matrix.csv`.

### Stage 2 - Distance & PLINK Prep (`step06`-`step10`)
- Convert IBS to dissimilarity, save genotype labels.
- Build PLINK `.map`, write TPED/TFAM, run PLINK `--r2` with LD decay R script.
- Outputs stored in `outputs/step06_convert_to_dissimilarity/`, `step07_extract_labels/`, `step08_create_plink_map/`, `step09_csv_to_tped/`, `step10_ld_decay/` (plot plus 50 kb bin summary with mean/median r^2). The analysis expects unfiltered PLINK LD output and only considers distances <= 5 Mb.

### Stage 3 - Pruning & Derived Datasets (`step11`-`step15`)
- Manual PowerShell helper executes PLINK pruning and recoding.
- Extract LD-pruned subset (`step12`), perform PCoA (`step13`), export STRUCTURE diploid (`step14`), standardize TASSEL PED (`step15`).
- Artefacts drive downstream manual (STRUCTURE) and TASSEL workflows.

### Stage 4 - STRUCTURE Analysis (Optional Manual)
- Copy `ld_pruned_genotypes_structure.txt` to `outputs/structure/raw/`.
- Run STRUCTURE for K = 1-10 (=10 reps), summarise with `structureHarvester.py` (harvester_input - harvester_output), align with CLUMPAK, record summary.

### Stage 5 - TASSEL GWAS (Manual GUI)
- Load `pruned_panel_fixed.ped` + `pruned_panel.map`.
- Export diagnostics (MAF, PCs, kinship) into `outputs/tassel/diagnostics/`.
- Import phenotypes, join covariates, run MLM/GLM (use Bonferroni threshold guidance).
- Save MLM stats (`outputs/tassel/mlm_stats.txt`), Manhattan/QQ plots (`outputs/tassel/plots/`), and document settings.

### Stage 6 - Post-GWAS Automation (`step16`-`step17`)
- `step16_trait_split`: split TASSEL stats by trait via R.
- `step17_trait_adjust`: apply Bonferroni and FDR corrections, produce top-N tables and `significant_snps_all_traits.csv`.
- Triggers overall summary (if configured) and feeds candidate gene review.

### Stage 7 - Candidate Gene Review & Reporting
- Merge adjusted trait outputs with probe map metadata for positional context.
- Cross-reference literature, update `outputs/reports/candidate_genes.xlsx`, and prioritise validation.

---

## 8. Outputs Checklist

1. `outputs/step03_filter_individuals/filtered_genotypes_strict.csv`
2. `outputs/step05_calculate_ibs/similarity_matrix.csv`
3. `outputs/step09_csv_to_tped/filtered_genotypes_strict.tped/.tfam`
4. `outputs/step10_ld_decay/ld_decay.png`, `.ld.summary`, `ld_decay_threshold.txt`, `ld_bins.tsv`, `ld_smooth.tsv`, `ld_decay_metrics.tsv`
5. `outputs/step11_plink_pruning/pruned_panel.{bed,bim,fam,ped,map}` + `.prune.in`
6. `outputs/step12_pruned_subset/ld_pruned_genotypes_strict.csv`
7. `outputs/step13_pcoa/coordinates.csv` + plot (`pcoa_plot.png`; pruned reruns write `pcoa_coordinates_pruned.csv` / `pcoa_plot_pruned.png`)
8. `outputs/step14_convert_to_structure/ld_pruned_genotypes_structure.txt`
9. `outputs/step15_fix_pruned_ped/pruned_panel_fixed.ped`
10. `outputs/tassel/mlm_stats.txt` + diagnostics/plots
11. `outputs/step17_trait_adjust/significant_snps_all_traits.csv`
12. `outputs/run_history/` log & progress JSON

---

## 9. Operational Notes

- **No per-step README files**: verify completion via console logs and generated artefacts.
- **Dry-run best practice**: run before executing subsets of steps to confirm configuration changes.
- **Manual stages** require consistent naming and output placement to keep automated stages functional.
- **Version control**: major changes to scripts or configuration must be reflected in both README and this overview to maintain alignment.

---

## 10. Quick Reference Commands

| Task | Command |
| --- | --- |
| Dry run (default config) | `python -m gwas_pipeline --dry-run` |
| Execute specific steps | `python -m gwas_pipeline --steps step01_filter_variants step02_filter_mapped_snps` |
| Run PLINK helper | `powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1` |
| Split TASSEL stats | `python -m gwas_pipeline --steps step16_trait_split` |
| Adjust traits | `python -m gwas_pipeline --steps step17_trait_adjust` |

---

## 11. Change Log Snapshot

- **Removed auto-generated step READMEs**: rely on logs/artefacts.
- **Structured STRUCTURE workspace**: `harvester_input/`, `harvester_output/`, `structureHarvester_scripts/`, `clumpak/`.
- **Added presentation flowchart** (`docs/flowchart.md`) mirroring README content.
- **Introduced post-GWAS candidate gene guidance and Bonferroni threshold documentation in README.**

---

This document is maintained alongside `README.md` and `.gitignore_context/copilot_instructions.md`. Updates to pipeline behaviour, directories, or tooling should be reflected here to preserve an accurate high-level reference for teammates and stakeholders.***



