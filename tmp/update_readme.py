from pathlib import Path

readme = Path("README.md")
text = readme.read_text(encoding="utf-8")
start = text.index("## Stage 4")
end = text.index("## Troubleshooting")
new = """## Stage 4 - PCoA and STRUCTURE input

**Automated pipeline steps**

- Run principal coordinates analysis on the pruned dissimilarity matrix (step13)
  `python -m gwas_pipeline --steps step13_pcoa`
  Uses `scripts/r/step13_pcoa.R` and writes `outputs/step13_pcoa/pcoa_coordinates_pruned.csv` plus `outputs/step13_pcoa/pcoa_plot_pruned.png`.

- Convert the LD-pruned genotype matrix to STRUCTURE diploid format (step14)
  `python -m gwas_pipeline --steps step14_convert_to_structure`
  Outputs `outputs/step14_convert_to_structure/ld_pruned_genotypes_structure.txt` (two-row diploid encoding per individual).

**Manual STRUCTURE runs (per SU-PBL spec)**

1. Copy `ld_pruned_genotypes_structure.txt` into your STRUCTURE workspace (recommended: `outputs/structure/raw/`).
2. Launch STRUCTURE (GUI or command-line) with the following baseline settings:
   - Admixture model ON, correlated allele frequencies.
   - Burn-in: 10,000 iterations; MCMC reps after burn-in: 10,000.
   - Test K values 1-10 with 10 replicates per K (set distinct random seeds).
   - Treat data as diploid (ploidy = 2) with `-9` as the missing genotype code.
3. Store each run under `outputs/structure/runs/INEW_<date>_K{value}_rep{n}/` to keep the run history alongside pipeline outputs.
4. Summarise runs with STRUCTURE HARVESTER (install separately) to compute Evanno delta K:
   ```powershell
   python structureHarvester_scripts\structureHarvester.py `
     --dir outputs\structure\runs `
     --out outputs\structure\harvester `
     --evanno
   ```
5. Zip the replicate folders per K (for example `K1.zip` ... `K10.zip`), combine them into `structure_results.zip`, and submit to https://clumpak.tau.ac.il/index.html for cluster alignment and barplots. Save the downloaded summaries back to `outputs/structure/clumpak/`.

Keep the HARVESTER `evanno.txt`/`summary.txt`, CLUMPAK consensus plots, and STRUCTURE logs with a short `_README.txt` noting the chosen K.

## Stage 5 - TASSEL GWAS preparation and execution

**Automated pipeline step**

- Normalise the PLINK PED for TASSEL (step15)
  `python -m gwas_pipeline --steps step15_fix_pruned_ped`
  Reads `outputs/step11_plink_pruning/pruned_panel.ped` and writes `outputs/step15_fix_pruned_ped/pruned_panel_fixed.ped` (sex column incremented by 1, allele `B` recoded to `C`).

**Prepare TASSEL inputs (PLINK)**

```powershell
plink --bfile outputs/step11_plink_pruning/pruned_panel `
      --extract outputs/step11_plink_pruning/pruned_panel.prune.in `
      --make-bed --allow-extra-chr `
      --out outputs/step11_plink_pruning/pruned_panel

plink --bfile outputs/step11_plink_pruning/pruned_panel `
      --recode `
      --allow-extra-chr `
      --tab `
      --out outputs/step11_plink_pruning/pruned_panel
```

This produces `outputs/step11_plink_pruning/pruned_panel.{bed,bim,fam,ped,map}` for TASSEL along with the pruned SNP list. If PLINK is not on PATH, call the executable explicitly (for example `& $env:PLINK_PATH --bfile ...`) or run the helper again to regenerate the panel.

**Recommended TASSEL MLM workflow (GUI)**

1. File -> Open As -> PLINK PED -> load `pruned_panel.ped`/`.map`.
2. Sequence View -> Data -> Geno Summary to inspect MAF, missingness, and heterozygosity.
3. Sequence View -> Analysis -> Relatedness -> Distance Matrix, then Analysis -> Relatedness -> MDS to generate principal coordinates (export eigenvalues and PCs).
4. Import phenotypes (tab-delimited) and any covariates that match sample IDs in the PED file.
5. Run the MLM (Mixed Linear Model) with:
   - Fixed effects: principal coordinates (PC1-PC3 from PCoA/MDS).
   - Kinship: TASSEL centered IBS or an external kinship matrix.
   - Evaluate GLM as a secondary model if MLM fails to converge.
6. Export per-trait MLM output via Analysis -> MLM -> Display Results -> Export. Save as `outputs/tassel/mlm_stats.txt`.

Archive Manhattan/Q-Q plots and TASSEL project files under `outputs/tassel/` with notes describing model terms.

## Stage 6 - Post-GWAS processing and reporting

- Split TASSEL MLM output per trait (step16):
  `python -m gwas_pipeline --steps step16_trait_split`
  Generates `outputs/step16_trait_split/traits/mlm_stats_<trait>.csv`.

- Apply multiple-testing correction and summarise top loci (step17):
  `python -m gwas_pipeline --steps step17_trait_adjust`
  Uses the `n_tests` (default 3048) and `top_n` (default 10) parameters from `config/pipeline.yaml` to create adjusted trait tables plus `outputs/step17_trait_adjust/significant_snps_all_traits.csv`.

After stage 6 completes the pipeline refreshes `outputs/summary/report.md` with key artifact paths; review and update the narrative (for example, chosen K and notable traits) before sharing results.
"""
updated = text[:start] + new + text[end:]
readme.write_text(updated, encoding="utf-8")
