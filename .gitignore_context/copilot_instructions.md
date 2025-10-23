# SU-PBL GWAS Pipeline — Copilot Instructions

These notes help Copilot (and other AI assistants) work safely within this repository. Keep the guidance aligned with the main `README.md` so human contributors and automation stay in sync.

---

## 1. Project Snapshot

- Hybrid Python/R/PLINK pipeline that progresses from raw SNP genotypes through QC, LD pruning, TASSEL GWAS, and post-GWAS reporting.
- Main orchestration entry point: `python -m gwas_pipeline`.
- Configuration: `config/pipeline.yaml` (primary run) and derivative configs for pruned workflows or demos.
- External executables required: PLINK (`plink.exe`) and Rscript.
- Per-step `_README.txt` files are **no longer generated automatically**; rely on logs and produced artefacts instead.

---

## 2. Repository Layout (unchanged but critical)

```
config/                  YAML files; keep structure stable.
data/                    Raw inputs and metadata (read-only for automation).
docs/                    Human documentation (Word spec, tutorials, demo data).
gwas_pipeline/           Thin wrapper so `python -m gwas_pipeline` works.
workflow/gwas_pipeline/  Core implementation (pipeline.py, steps.py, etc.).
scripts/python/          Utility scripts (self-check, validation helpers).
scripts/r/               R scripts invoked by pipeline steps.
scripts/cli_examples/    PowerShell helpers (e.g., PLINK pruning wrapper).
outputs/                 Step outputs (`outputs/stepNN_*`), run logs in `outputs/run_history/`.
outputs/structure/       STRUCTURE workspace (`raw/`, `runs/`, `structure_harvester/{harvester_input,harvester_output,structureHarvester_scripts}`, `clumpak/`).
```

Maintain path casing and naming conventions; downstream tooling expects them verbatim.

---

## 3. Key Execution Facts

- `workflow/gwas_pipeline/pipeline.py` defines the canonical step order (`STEP_SEQUENCE`) and handles run logging.
- `workflow/gwas_pipeline/steps.py` implements each step. Many steps shell out to PLINK or R. Respect existing helper functions (e.g., `step_logger`, `_run_r_script`).
- Since `_write_step_summary` and `_generate_overall_summary` were removed, avoid reintroducing auto-generated per-step README files unless the team agrees to reinstate them.
- Logs: every run creates `outputs/run_history/run_YYYYMMDD_HHMMSS{_dryrun}.log/json`. Use these for diagnostics instead of relying on text files in each step directory.
- TASSEL work (manual GUI) must write outputs back under `outputs/tassel/` so Steps 16–17 can continue automatically.

---

## 4. Implementation Guardrails

1. **Do not rename steps** – step identifiers (`step01_filter_variants`, etc.) appear in configs, scripts, and documentation.
2. **Read configuration values** – hard-code nothing that already lives in `pipeline.yaml`. Leverage `PipelineConfig.path()` helpers.
3. **Validate external tools gracefully** – when adding features that call PLINK/R, copy the defensive patterns already used in Step 10.
4. **Keep docs in sync** – any change to behaviour or required paths must be reflected both here and in `README.md`.
5. **Preserve Windows compatibility** – PowerShell commands and path handling must work on Windows first, while remaining portable where possible.
6. **No silent deletions** – if you retire a file or folder, update documentation and scripts that reference it.

---

## 5. Stage Overview (for quick reference)

| Stage | Steps | Summary | Notes |
|---|---|---|---|
| Stage 1 | 01-05 | SNP-level QC, individual filtering, matrix transpose, IBS similarity. | Produces strict genotype CSV and similarity matrix. |
| Stage 2 | 06-10 | Dissimilarity, label extraction, PLINK map/TPED conversion, LD decay analysis. | Requires PLINK and Rscript. Step 10 writes the plot, bins, and metrics to `outputs/step10_ld_decay/`.
| Stage 3 | 11-15 | Manual PLINK pruning helper, pruned subset export, PCoA (unpruned + pruned reruns), STRUCTURE conversion, TASSEL-ready PED fix. | Step 11 manual: `scripts/cli_examples/run_plink_stage.ps1`; STRUCTURE input created in Step 14. |
| Stage 4 | STRUCTURE (optional) | Manual runs, HARVESTER summaries, CLUMPAK alignment. | Use `outputs/structure/structure_harvester/{harvester_input,harvester_output,structureHarvester_scripts}` and `outputs/structure/clumpak/`. |
| Stage 5 | TASSEL GUI | Import genotypes/phenotypes, run diagnostic exports, compute MLM. | Outputs must be saved under `outputs/tassel/`; remind users about Bonferroni thresholds. |
| Stage 6 | 16-17 | Trait split and multiple-testing adjustment via R scripts. | Consumes TASSEL outputs; generates per-trait and cross-trait summaries. |

When adding new functionality, extend this table and the README simultaneously.

---

## 6. Manual Workflow Reminders

- The pipeline no longer writes `_README.txt` per step. If human-readable notes are required, authors must add them manually (e.g., TASSEL log).
- Encourage users to consult the latest log (`outputs/run_history/latest.json`) to confirm completion status.
- `helper_self_check.py` is the preferred pre-flight script; mention it when guiding users through setup.
- STRUCTURE tooling now lives under `outputs/structure/structure_harvester/{harvester_input,harvester_output,structureHarvester_scripts}` with CLUMPAK summaries stored alongside in `outputs/structure/clumpak/`; keep instructions aligned with that layout.

---

## 7. Common Pitfalls & Messaging

- **Missing executables**: instruct users to set `$env:PLINK_PATH` or extend PATH with the R `bin` directory before re-running Step 10.
- **Incorrect outputs directory**: remind users that steps assume `outputs/stepNN_<name>/`; altering these paths requires YAML updates.
- **TASSEL exports misplaced**: emphasise saving results back into the repository tree, otherwise Steps 16–17 fail.
- **Non-technical audience**: when generating helper text or prompts, prefer plain language and explain why each command runs.

---

## 8. Testing & Verification Guidance

- Python code: `python -m py_compile workflow/gwas_pipeline/steps.py` is a lightweight syntax check.
- R scripts: ensure they stay compatible with command-line execution (`Rscript script.R ...`).
- When instructing users, recommend running `python -m gwas_pipeline --dry-run` after configuration changes to catch missing paths early.

---

## 9. When to Escalate

Contact maintainers (or leave comments) if:

- Introducing new stages or altering step order.
- Modifying how logs are written (could break downstream monitoring).
- Adding dependencies beyond what `workflow/requirements.txt` captures.
- Reverting to automatic README generation – requires consensus and documentation updates.

Stay disciplined about mirroring behavioural changes in both the README and these instructions. This keeps newcomers and automation aligned with the actual behaviour of the pipeline.
