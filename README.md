# SU-PBL GWAS Pipeline

A reproducible, stepwise workflow that takes raw SNP microarray genotypes through quality control, LD-based pruning, population-structure analyses, GWAS association testing, and post-GWAS reporting. The numbered steps are configured by `config/pipeline.yaml` and orchestrated by `python -m gwas_pipeline`.

---

## 0. Workspace Setup (do this once per machine)

| Step | Command (PowerShell) | Why it matters |
|---|---|---|
| 1. Create an isolated Python environment | `python -m venv .venv` <br> `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass` <br> `.\\.venv\\Scripts\\Activate.ps1` | Keeps dependencies for the pipeline separate from global Python installs. The execution-policy change is temporary and permits the activation script to run; activate this environment before each pipeline session. |
| 2. Install Python requirements | `pip install -r workflow/requirements.txt` | Installs `pandas`, `numpy`, `PyYAML`, and supporting libraries used by the Python steps. |
| 3. Verify external executables | Confirm `plink` and `Rscript` run from the terminal (`plink --version`, `Rscript --version`). If not, set:<br>`$env:PLINK_PATH = "C:\\Path\\to\\plink.exe"`<br>`$env:Path = "C:\\Program Files\\R\\R-4.5.0\\bin;$env:Path"` | Several steps shell out to PLINK and R. Setting the environment variables avoids configuration errors later. |
| 4. Copy or link input data | Place genotype CSVs and supporting metadata in `data/`. Sample files are available in `docs/sample_data/`. | Ensures the relative paths used in `pipeline.yaml` are valid. |
| 5. Sanity check | `python scripts/python/helper_self_check.py` | Confirms required inputs exist and highlights anything missing before you attempt a run. |

Keep PowerShell open with the virtual environment activated while you work (`.\\.venv\\Scripts\\Activate.ps1`). Deactivate with `deactivate` when finished.

---

## 1. Directory Overview

| Folder | Contents & guidance |
|---|---|
| `config/` | YAML configuration files. `pipeline.yaml` drives the main workflow; `pipeline.sample.yaml` is bundled for practice. Keep file names unchanged so automation finds them. |
| `data/raw/` | Primary genotype matrix (CSV) plus any phenotype tables used later. Treat as read-only once validated. |
| `data/metadata/` | Probe annotation resources such as the 35k Excel sheet (`35k_probe_set_IWGSCv1.xlsx`) and the cleaned map (`35k_probe_set_IWGSCv1_cleaned_tab.map`). |
| `docs/` | Human-facing documentation including the Word specification and tutorial materials. Demo inputs live under `docs/sample_data/`. |
| `gwas_pipeline/` | Lightweight module so you can run `python -m gwas_pipeline`. Do not modify unless you are changing the CLI interface. |
| `workflow/gwas_pipeline/` | Core Python implementation (`pipeline.py`, `steps.py`, `config.py`, `progress.py`, `utils.py`). All automated steps originate here. |
| `scripts/python/` | Helper utilities (`helper_self_check.py`, validation tools). Extend these for diagnostics rather than editing the core pipeline. |
| `scripts/r/` | R scripts executed by certain steps (`step10_ld_decay_enhanced.R`, `step13_pcoa.R`, etc.). |
| `scripts/cli_examples/` | PowerShell helpers, e.g., `run_plink_stage.ps1` which wraps the PLINK pruning sequence. |
| `outputs/` | The pipeline writes each step to `outputs/stepNN_<name>/`. Logs and history land in `outputs/run_history/`. No additional step-level README files are generated�inspect the files themselves to verify completion. |
| `outputs/structure/` | STRUCTURE workspace organised as: raw exports in `raw/`, each STRUCTURE run under `runs/<DATE>_K#/`, HARVESTER inputs in `structure_harvester/harvester_input/`, summaries in `structure_harvester/harvester_output/`, helper scripts in `structure_harvester/structureHarvester_scripts/`, and CLUMPAK bundles under `clumpak/` (`clumpak/input/` for zipped runs, `clumpak/outputs/` for consensus plots). |
| `outputs/tassel/` | TASSEL workspace: `diagnostics/` (with `genotype_summary/`, `kinship/`, `mds/`, `mlm_input_table.txt`), `phenotypes/`, `plots/`, `derived/`, `project/`, top-level `mlm_stats.txt`, `mlm_effects.txt`, and an optional manually maintained `_README.txt`. |

---

## 2. Key Inputs and Expected Formats

| Item | Path (relative) | Description & notes |
|---|---|---|
| Genotype matrix | `data/raw/genotypes/raw_genotypes.csv` | Rows = SNP probes, columns = `Probe ID` + individual genotype codes (0/1/2/9). |
| Probe annotation (Excel) | `data/metadata/35k_probe_set_IWGSCv1.xlsx` | Used to retain SNPs with mapped chromosomes. Must contain columns `35K SNPId` and `IWGSC_v1_Chromosome`. |
| Probe annotation (map) | `data/metadata/35k_probe_set_IWGSCv1_cleaned_tab.map` | Tab-delimited file with cleaned coordinates for PLINK `.map` creation. |
| Phenotype table (for TASSEL) | `outputs/tassel/phenotypes/<run_label>.txt` (you create this) | First column `Taxa`, remaining columns = traits. TASSEL consumes this during MLM. |
| Optional configs | `config/pipeline.pruned.yaml` etc. | Provide alternative output folders or parameters for pruned workflows. |

Check `config/pipeline.yaml` to confirm these paths match your data layout. Most step blocks already reference these defaults. Adjust with care and re-run `helper_self_check.py` after edits.

---

## 3. Running a Sample Workflow (safe sandbox)

The repository includes a miniature dataset so you can practice without touching real project data.

```powershell
# List available steps and confirm the configuration loads
python -m gwas_pipeline --config config/pipeline.sample.yaml --dry-run

# Execute the full sample flow (writes to outputs/sample_run/stepNN_*)
python -m gwas_pipeline --config config/pipeline.sample.yaml
```

Inspect `outputs/sample_run/` to understand the artifacts created at each step, then apply the same commands (without the `--config` override) to your real dataset once you are comfortable.

---

## 4. Pipeline Stages (main run)

The table below summarises each numbered step, the command to run it on its own, required inputs, outputs, and the purpose. Unless stated, run commands from the project root with your virtual environment active.

### Stage 1 - Microarray QC and Formatting (Steps 01-05)

| Step | Command | Key inputs | Outputs | Purpose & context |
|---|---|---|---|---|
| 01 `filter_variants` | `python -m gwas_pipeline --steps step01_filter_variants` | `data/raw/genotypes/raw_genotypes.csv` | `outputs/step01_filter_variants/filtered_genotypes.csv`, removal report | Removes SNPs failing missingness (default 5%) or minor allele frequency thresholds to stabilise downstream analyses. |
| 02 `filter_mapped_snps` | `python -m gwas_pipeline --steps step02_filter_mapped_snps` | Step 01 output + `35k_probe_set_IWGSCv1.xlsx` | `outputs/step02_filter_mapped_snps/filtered_genotypes_mapped.csv`, SNP keep list | Retains only probes with chromosome assignments so physical positions exist for PLINK and LD decay. |
| 03 `filter_individuals` | `python -m gwas_pipeline --steps step03_filter_individuals` | Step 02 output | `outputs/step03_filter_individuals/filtered_genotypes_strict.csv`, removed individuals CSV | Removes samples exceeding the missing-call threshold (default 5%) to prevent noisy genotype vectors. |
| 04 `transpose_matrix` | `python -m gwas_pipeline --steps step04_transpose_matrix` | Step 03 output | `outputs/step04_transpose_matrix/filtered_genotypes_strict_T.csv` | Rotates SNP rows into columns so each row represents an individualâ€"this orientation is required by PLINK conversion and similarity calculations. |
| 05 `calculate_ibs` | `python -m gwas_pipeline --steps step05_calculate_ibs` | Step 04 output | `outputs/step05_calculate_ibs/similarity_matrix.csv` | Computes the Identity-by-State similarity matrix between individuals, handling missing data gracefully. |

**Checkpoint:** verify the strict genotype CSV and similarity matrix exist before moving on (`outputs/step03_filter_individuals/`, `outputs/step05_calculate_ibs/`). The logs live in `outputs/run_history/` if you need timing or error details.

### Stage 2 - Distance Matrices and PLINK Preparation (Steps 06-10)

| Step | Command | Key inputs | Outputs | Purpose & context |
|---|---|---|---|---|
| 06 `convert_to_dissimilarity` | `python -m gwas_pipeline --steps step06_convert_to_dissimilarity` | Step 05 similarity matrix | `outputs/step06_convert_to_dissimilarity/dissimilarity_matrix.csv` | Converts similarity to distance (1 âˆ’ IBS), enabling ordination and clustering in later stages. |
| 07 `extract_labels` | `python -m gwas_pipeline --steps step07_extract_labels` | Step 05 similarity matrix | `outputs/step07_extract_labels/genotype_labels.txt` | Captures the ordered list of genotype IDs for downstream R scripts (PCoA). |
| 08 `create_plink_map` | `python -m gwas_pipeline --steps step08_create_plink_map` | SNP keep list + cleaned map file | `outputs/step08_create_plink_map/filtered_genotypes_strict.map` | Builds a PLINK `.map` with chromosome and base-pair positions so PLINK understands genome coordinates. |
| 09 `csv_to_tped` | `python -m gwas_pipeline --steps step09_csv_to_tped` | Transposed genotypes + PLINK map | `outputs/step09_csv_to_tped/filtered_genotypes_strict.tped` and `.tfam` | Translates the genotype matrix into PLINK TPED/TFAM format, performing 0/1/2/9 to allele-pair conversion. |

Before running Step 10, confirm PLINK and Rscript are available to the current PowerShell session (skip if already on PATH):

```powershell
# Set PLINK to the full installed executable on this machine
$env:PLINK_PATH = "C:\Program Files\PLINK\plink_win64_20241022\plink.exe"
# Optionally add R to PATH for this session
$env:Path = "C:\Program Files\R\R-4.5.0\bin;$env:Path"
```

Use `plink --version` and `Rscript --version` to verify the commands resolve.

| Step | Command | Key inputs | Outputs | Purpose & context |
|---|---|---|---|---|
| 10 `ld_decay_plot` | `python -m gwas_pipeline --steps step10_ld_decay_plot` | TPED/TFAM + `scripts/r/step10_ld_decay_enhanced.R` | `outputs/step10_ld_decay/filtered_genotypes_strict.ld.summary`, `.ld`, `ld_decay.png`, `ld_decay_threshold.txt`, `ld_bins.tsv`, `ld_smooth.tsv`, `ld_decay_metrics.tsv` | Runs `plink --r2`, reduces the LD output to distance vs r^2, then bins the distances in 50 kb increments to report mean/median r^2 and plot LD decay. |

**Before continuing:** check that the LD decay plot renders without errors and review `ld_bins.tsv`, `ld_smooth.tsv`, and `ld_decay_metrics.tsv` to confirm the mean/median r^2 trend across distance bins. The compatibility file `ld_decay_threshold.txt` will show `NA` values when no explicit r^2 crossing is calculated. 

### Stage 3 - LD Pruning and Pruned Dataset (Steps 11-15)

| Step | How to run | Inputs expected | Outputs | Purpose & context |
|---|---|---|---|---|
| 11 `plink_pruning` | Manual helper: `powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1` (run after Step 10) | TPED/TFAM from Step 09, `.map` from Step 08 | `outputs/step11_plink_pruning/pruned_panel.*`, `.prune.in`, log files | Executes a sequence of PLINK commands: LD pruning (`--indep-pairwise`), QC binary conversion, pruned binary export, and PED/MAP generation for TASSEL. Review PowerShell output to confirm success. |
| 12 `extract_pruned_subset` | `python -m gwas_pipeline --steps step12_extract_pruned_subset` | `.prune.in` list + strict genotype CSV | `outputs/step12_pruned_subset/ld_pruned_genotypes_strict.csv` | Filters the original genotype matrix to only the pruned SNPs, producing a manageable CSV for STRUCTURE or other analyses. |
| 13 `pcoa` | `python -m gwas_pipeline --steps step13_pcoa` | Dissimilarity matrix (Step 06) + genotype labels (Step 07) + `scripts/r/step13_pcoa.R` | `outputs/step13_pcoa/coordinates.csv`, `pcoa_plot.png` | Performs Principal Coordinates Analysis via R to visualise population structure and identify clustering. |
| 14 `convert_to_structure` | `python -m gwas_pipeline --steps step14_convert_to_structure` | Pruned genotype CSV (Step 12) | `outputs/step14_convert_to_structure/ld_pruned_genotypes_structure.txt` | Converts the 0/1/2 coding into STRUCTURE's two-column allele format; ready for STRUCTURE or similar tools. |
| 15 `fix_pruned_ped` | `python -m gwas_pipeline --steps step15_fix_pruned_ped` | Pruned PED from Step 11 | `outputs/step15_fix_pruned_ped/pruned_panel_fixed.ped` | Normalises PED formatting for TASSEL (sex codes, allele naming). This file pairs with `pruned_panel.map` for TASSEL import. |

If you maintain a pruned-only configuration (`config/pipeline.pruned.yaml`), you can re-run Steps 04-07 and 13 with the pruned inputs using the suffixed step names (`step04_transpose_pruned`, etc.). This keeps results separated from the main outputs.

```powershell
# Refresh the baseline (unpruned) matrices and PCoA
python -m gwas_pipeline --steps `
  step04_transpose_matrix `
  step05_calculate_ibs `
  step06_convert_to_dissimilarity `
  step07_extract_labels `
  step13_pcoa
```
(Writes `outputs/step13_pcoa/coordinates.csv` and `pcoa_plot.png`.)

```powershell
# Launch the LD pruning helper (produces pruned_panel.* plus prune lists)
powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1

# Example: regenerate the pruned-only artefacts
python -m gwas_pipeline --config config/pipeline.pruned.yaml --steps `
  step04_transpose_pruned `
  step05_calculate_ibs_pruned `
  step06_convert_to_dissimilarity_pruned `
  step07_extract_labels_pruned `
  step13_pcoa
```
(Run the pruning helper above once Step 10 has produced the TPED/TFAM inputs.)
(Produces `outputs/step13_pcoa/pcoa_coordinates_pruned.csv` and `pcoa_plot_pruned.png` without touching the unpruned files.)

### Stage 4 - STRUCTURE Analysis (Optional Manual Workflow)

Use the STRUCTURE-ready file from Step 14 if you need population clustering beyond PCoA:

1. Copy `outputs/step14_convert_to_structure/ld_pruned_genotypes_structure.txt` into your STRUCTURE workspace (recommended destination: `outputs/structure/raw/`).
2. Launch STRUCTURE (GUI or command line) with these baseline settings:
   - Admixture model **ON**, correlated allele frequencies.
   - Burn-in: 10000 iterations; MCMC reps after burn-in: 10000.
   - Test K values 1-10 with 10 replicates per K (set distinct random seeds).
   - Treat data as diploid (`ploidy = 2`) with `-9` as the missing genotype code.
3. Store each run under `outputs/structure/runs/<DATE>_K<value>_rep<n>/` so the history sits alongside the pipeline outputs. Copy or symlink the runs you wish to summarise into `outputs/structure/structure_harvester/harvester_input/`.
4. (Optional but recommended) Summarise runs with STRUCTURE HARVESTER (stored under `outputs/structure/structure_harvester/structureHarvester_scripts/`) to compute Evanno Delta K:

   ```powershell
   python outputs/structure/structure_harvester/structureHarvester_scripts/structureHarvester.py `
     --dir outputs/structure/structure_harvester/harvester_input `
     --out outputs/structure/structure_harvester/harvester_output `
     --evanno
   ```

5. Zip the replicate folders per K (for example, `K1.zip to K10.zip`), combine them into `structure_results.zip`, and submit to [CLUMPAK](https://clumpak.tau.ac.il/index.html) for cluster alignment and barplots. Store the downloads in `outputs/structure/clumpak/`.
6. Retain HARVESTER files (`evanno.txt`, `summary.txt`) inside `outputs/structure/structure_harvester/harvester_output/`, keep CLUMPAK consensus plots under `outputs/structure/clumpak/`, and add a short `_README.txt` noting the chosen K and key observations.

### Stage 5 - TASSEL-Based GWAS (Manual GUI Workflow)

TASSEL work happens outside Python but should obey the same directory conventions so automated post-processing can find the files.

1. **Load genotypes**: `File -> Open As -> PLINK Ped`, selecting `outputs/step15_fix_pruned_ped/pruned_panel_fixed.ped` and `outputs/step11_plink_pruning/pruned_panel.map`. Optionally save the TASSEL project under `outputs/tassel/project/`.
2. **Diagnostic exports** (all under `outputs/tassel/diagnostics/`): Geno summary tables, MDS eigenvalues, MDS coordinates, and the centered IBS kinship matrix.
3. **Phenotype import**: store the source table under `outputs/tassel/phenotypes/<run_label>.txt`, then load it in TASSEL.
4. **Join covariates**: `Data -> Intersect Join` to combine phenotypes, PCs, and any covariates. Export the joined table as `mlm_input_table.tsv` in diagnostics.
5. **Run MLM (or GLM fallback)**: `Analysis -> Associations -> MLM`, using the kinship file and PCs as fixed effects. Export statistics to `outputs/tassel/mlm_stats.txt`. Save Manhattan and QQ plots to `outputs/tassel/plots/`.
6. **Record notes**: create or update `outputs/tassel/_README.txt` manually with TASSEL version, chosen covariates, and any issues encountered.

The key hand-off file for the automated post-GWAS stage is `outputs/tassel/mlm_stats.txt`.

**Bonferroni threshold reminder**  
The classical Bonferroni cut-off for TASSEL p-values is `p <= 1 / N` (a = 0.05), where `N` is the total number of markers tested. Expressed on the `-log10` scale used in Manhattan plots, the threshold becomes `-log10(0.05 / X)` with `X = N`. Use this value when interpreting TASSEL plots and when configuring downstream filtering.

### Stage 6 - Post-GWAS Processing (Steps 16-17)

| Step | Command | Inputs | Outputs | Purpose & context |
|---|---|---|---|---|
| 16 `trait_split` | `python -m gwas_pipeline --steps step16_trait_split` | TASSEL MLM stats (`outputs/tassel/mlm_stats.txt`) + `scripts/r/step16_trait_split.R` | `outputs/step16_trait_split/traits/mlm_stats_<trait>.csv` | Splits the TASSEL results into one CSV per trait so adjustments can be applied independently. |
| 17 `trait_adjust` | `python -m gwas_pipeline --steps step17_trait_adjust` | Trait CSV directory from Step 16 + `scripts/r/step17_trait_adjust.R` | `outputs/step17_trait_adjust/traits/adj_p_<trait>.csv`, `top_<n>_<trait>.csv`, `significant_snps_all_traits.csv` | Applies Bonferroni and FDR corrections, highlights the top SNPs per trait, and collates significant findings into a single summary. |

After Step 17, review the generated tables and plots, document biological interpretations in your lab notes, and prepare any reports needed for collaborators.

### Candidate Gene Identification and Reporting

Use the adjusted trait outputs to shortlist candidate loci and connect them to biological context:

1. **Start with summary tables.** Open `outputs/step17_trait_adjust/significant_snps_all_traits.csv` to review every SNP surpassing the multiple-testing threshold (Bonferroni by default). For trait-level inspection, pair each `adj_p_<trait>.csv` (full statistics) with its `top_<n>_<trait>.csv` (prioritised subset).
   - At this stage you are defining credible hits. Remember that LD can cluster signals—multiple significant SNPs may tag the same locus. Perform LD clumping or identify independent signals before final ranking.
2. **Cross-reference functional resources.** Using the SNP coordinates, query genome/annotation databases (IWGSC, EnsemblPlants, or crop-specific browsers). Capture known trait associations, functional annotations (GO, protein domains), and expression context.
   - At a minimum annotate each locus with SNP ID, chromosome/position, nearest gene(s), gene ID/function, SNP–gene distance, and any prior evidence (expression, literature, regulatory data).
3. **Summarise findings.** Track candidates and supporting evidence in `outputs/reports/candidate_genes.xlsx` (or your team’s tracker). Link back to the trait CSVs, include any regional association plots (e.g., TASSEL output) or gene model snapshots, and note the rationale for prioritisation.
   - Recommended columns: Trait, SNP, p-value, effect size, chromosome/position, implicated gene(s), annotation, prior evidence, proposed follow-up action.
4. **Plan validation.** For each locus, record whether replication is pending, marker development (e.g., KASP/PCR) is underway, or functional validation (expression assays, genome editing, field trials) is scheduled.
   - Maintain a “validation status” entry (e.g., “Pending replication”, “Marker design in progress”, “Functional assay planned 2026”) so downstream teams can prioritise work.

Keeping these notes and files in-repo ensures collaborators understand how each candidate was selected, what evidence supports the decision, and which artefacts underpin prioritisation.

---

## 5. Troubleshooting & Validation

| Issue | Likely cause | Resolution |
|---|---|---|
| `ConfigError: Missing input path` | Path not populated in `pipeline.yaml` | Update the relevant `steps.<name>` block or `paths` entry, then rerun. |
| `PLINK executable not found` | PLINK not on PATH | Set `$env:PLINK_PATH` or configure `tools.plink_executable` in `pipeline.yaml`. |
| `Rscript not found` | R not installed or PATH missing | Install R, then update PATH or `tools.rscript_executable`. |
| Step 10 fails with empty LD results | Dataset too small or map coordinates missing | Recheck `.map` file, confirm Step 08 succeeded, and verify PLINK log for warnings. |
| TASSEL outputs missing | Files exported to different folders | Re-export into the expected `outputs/tassel/` subdirectories so Steps 16-17 can locate them. |
| Want to verify run artifacts | Need quick checklist | Use the list below or run an internal validation script if available. |

### Quick Verification Checklist

Ensure the following files exist before engaging in interpretation or sharing results:

- `outputs/step03_filter_individuals/filtered_genotypes_strict.csv`
- `outputs/step08_create_plink_map/filtered_genotypes_strict.map`
- `outputs/step09_csv_to_tped/filtered_genotypes_strict.{tped,tfam}`
- `outputs/step10_ld_decay/ld_decay.png`, `.ld.summary`, `ld_decay_threshold.txt`, `ld_bins.tsv`, `ld_smooth.tsv`, `ld_decay_metrics.tsv`
- `outputs/step11_plink_pruning/pruned_panel.{bed,bim,fam,ped,map}` and `.prune.in`
- `outputs/step12_pruned_subset/ld_pruned_genotypes_strict.csv`
- `outputs/step13_pcoa/pcoa_plot.png`
- `outputs/step14_convert_to_structure/ld_pruned_genotypes_structure.txt`
- `outputs/step15_fix_pruned_ped/pruned_panel_fixed.ped`
- `outputs/tassel/mlm_stats.txt`
- `outputs/step16_trait_split/traits/` (per-trait CSVs)
- `outputs/step17_trait_adjust/significant_snps_all_traits.csv`

Logs for each execution sit in `outputs/run_history/`. Open the latest JSON or log file if you need timestamps, durations, or error messages.

---

## 6. Suggested Working Routine

1. **Activate environment**: `.\\.venv\\Scripts\\Activate.ps1`  
2. **Health check**: `python scripts/python/helper_self_check.py`  
3. **Dry run** (optional but recommended when editing configs): `python -m gwas_pipeline --dry-run`  
4. **Run stages**: either run everything (`python -m gwas_pipeline`) or call specific steps as needed.  
5. **Manual tasks**: execute Step 11 via the PowerShell helper, complete TASSEL workflow, and store exports in the agreed folders.  
6. **Post-processing**: run Steps 16-17, review outputs, and take notes.  
7. **Clean up**: deactivate the environment with `deactivate` when finished.  

Document major decisions (e.g., LD threshold choice, TASSEL covariates) in your lab notebook to maintain traceability.

---

## 7. Need Assistance<

If you encounter unexpected behaviour:

- Review the most recent log in `outputs/run_history/` for the exact error message.
- Re-run the individual step with `--dry-run` first to confirm configuration, then execute it normally.
- Use the sample pipeline to confirm toolchain installation (`python -m gwas_pipeline --config config/pipeline.sample.yaml --steps ...`). If the sample fails, prerequisite tools are likely misconfigured.
- Share relevant log snippets and configuration values when requesting help from teammates; this speeds up troubleshooting.