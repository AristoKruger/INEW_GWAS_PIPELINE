from __future__ import annotations

import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd

from .config import ConfigError, PipelineConfig
from .utils import ensure_parent, get_logger, step_logger, write_text


logger = get_logger()


def _stringify_missing_code(value: object) -> str:
    return str(value) if value is not None else "9"


def _resolve_input_path(config: PipelineConfig, step_cfg: dict, key: str, fallback_keys: Iterable[str] | None = None) -> Path:
    raw = step_cfg.get(key)
    if raw:
        path = config.resolve_path(raw)
        if path is None:
            raise ConfigError(f"Input path for {key} is not set")
        return path
    if fallback_keys:
        fallback = config.path(*fallback_keys)
        if fallback:
            return fallback
    raise ConfigError(f"Missing input path for key '{key}' and no fallback provided")


def _to_numeric_dataframe(df: pd.DataFrame, missing_code: str | int) -> pd.DataFrame:
    numeric = df.apply(pd.to_numeric, errors="coerce")
    if missing_code is not None:
        numeric = numeric.replace(float(missing_code), np.nan)
    return numeric


def _load_step_config(config: PipelineConfig, primary: str, legacy: str | None = None):
    keys = [primary]
    if legacy and legacy != primary:
        keys.append(legacy)
    for key in keys:
        step_cfg = config.get("steps", key)
        if isinstance(step_cfg, dict):
            return step_cfg, key
    return None, primary




def plink_pruning(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step11_plink_pruning", "plink_pruning")
    if step_cfg is None:
        logger.info("Skipping step11_plink_pruning (not configured)")
        return

    panel_prefix = config.path("steps", step_key, "panel_prefix", create_parent=True)
    if panel_prefix is None:
        raise ConfigError("steps.step11_plink_pruning.panel_prefix must be set")

    output_dir = panel_prefix.parent
    ensure_parent(output_dir / "_README.txt")
    message = "\n".join(
        [
            "Run PLINK LD pruning manually to generate the binary and PED/MAP panel.",
            "",
            "Suggested command (PowerShell):",
            "  powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1",
            "",
            f"Outputs should be written under: {output_dir}",
        ]
    )
    write_text(output_dir / "_README.txt", message)
    logger.info("Manual PLINK pruning required; see %s for instructions.", output_dir / "_README.txt")


def filter_variants(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step01_filter_variants", "filter_variants")
    if step_cfg is None:
        logger.info("Skipping step01_filter_variants (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("paths", "raw_genotypes"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step01_filter_variants.output_csv must be set")

    removed_log_path = config.path("steps", step_key, "removed_log", create_parent=True)
    maf_cutoff = float(step_cfg.get("maf_cutoff_percent", 5.0)) / 100.0
    missing_cutoff = float(step_cfg.get("missing_cutoff_percent", 10.0)) / 100.0
    id_column = step_cfg.get("snp_id_column", "Probe ID")
    missing_code = _stringify_missing_code(step_cfg.get("missing_code", 9))

    with step_logger("Filter SNPs by MAF and missingness"):
        df = pd.read_csv(input_path)
        if id_column not in df.columns:
            raise ConfigError(f"Column '{id_column}' not found in {input_path}")

        geno = df.drop(columns=[id_column])
        numeric = geno.apply(pd.to_numeric, errors="coerce")
        numeric = numeric.replace(float(missing_code), np.nan)

        missing_fraction = numeric.isna().mean(axis=1)
        called_genotypes = numeric.shape[1] - numeric.isna().sum(axis=1)
        called_alleles = called_genotypes * 2

        count_alt = (numeric == 0).sum(axis=1) * 2 + (numeric == 1).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            alt_freq = count_alt / called_alleles
        maf = np.minimum(alt_freq, 1 - alt_freq)
        maf = maf.fillna(0)

        keep_mask = (
            (called_alleles > 0)
            & (missing_fraction <= missing_cutoff)
            & (maf >= maf_cutoff)
        )

        filtered = df.loc[keep_mask]
        ensure_parent(output_path)
        filtered.to_csv(output_path, index=False)

        removed = df.loc[~keep_mask, [id_column]].copy()
        removed["missing_fraction"] = missing_fraction.loc[~keep_mask].round(4)
        removed["maf"] = maf.loc[~keep_mask].round(4)
        reasons: List[str] = []
        for idx in removed.index:
            badges: List[str] = []
            if called_alleles.loc[idx] == 0:
                badges.append("no_calls")
            if missing_fraction.loc[idx] > missing_cutoff:
                badges.append(f"missing>{missing_cutoff * 100:.2f}%")
            if maf.loc[idx] < maf_cutoff:
                badges.append(f"maf<{maf_cutoff * 100:.2f}%")
            reasons.append(";".join(badges) or "filtered")
        removed["reason"] = reasons
        if removed_log_path is not None and not removed.empty:
            ensure_parent(removed_log_path)
            removed.to_csv(removed_log_path, index=False)

        logger.info(
            "Filtered SNPs: kept %d of %d (MAF = %.2f%%, missing = %.2f%%)",
            keep_mask.sum(),
            len(df),
            maf_cutoff * 100,
            missing_cutoff * 100,
        )


def filter_mapped_snps(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step02_filter_mapped_snps", "filter_mapped_snps")
    if step_cfg is None:
        logger.info("Skipping step02_filter_mapped_snps (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step01_filter_variants", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step02_filter_mapped_snps.output_csv must be set")

    snp_list_path = config.path("steps", step_key, "snp_list_txt", create_parent=True)
    if snp_list_path is None:
        raise ConfigError("steps.step02_filter_mapped_snps.snp_list_txt must be set")

    annotation_path = _resolve_input_path(
        config,
        step_cfg,
        "probe_annotation_excel",
        fallback_keys=("paths", "probe_annotation_excel"),
    )

    id_column = step_cfg.get("snp_id_column", "Probe ID")
    snp_id_column = step_cfg.get("annotation_snp_column", "35K SNPId")
    chromosome_column = step_cfg.get("chromosome_column", "IWGSC_v1_Chromosome")

    with step_logger("Filter SNPs with IWGSC chromosome assignments"):
        geno_df = pd.read_csv(input_path)
        if id_column not in geno_df.columns:
            raise ConfigError(f"Column '{id_column}' not found in {input_path}")

        probe_df = pd.read_excel(annotation_path, dtype=str)
        valid_snps = probe_df.loc[probe_df[chromosome_column].notna(), snp_id_column].astype(str).unique()
        mask = geno_df[id_column].astype(str).isin(valid_snps)
        filtered = geno_df.loc[mask]
        ensure_parent(output_path)
        filtered.to_csv(output_path, index=False)

        ensure_parent(snp_list_path)
        filtered[id_column].to_csv(snp_list_path, index=False, header=False)

        logger.info(
            "Mapped SNPs: kept %d of %d using annotation %s",
            mask.sum(),
            len(geno_df),
            annotation_path,
        )


def filter_individuals(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step03_filter_individuals", "filter_individuals")
    if step_cfg is None:
        logger.info("Skipping step03_filter_individuals (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step02_filter_mapped_snps", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step03_filter_individuals.output_csv must be set")

    removed_log_path = config.path("steps", step_key, "removed_log", create_parent=True)
    id_column = step_cfg.get("snp_id_column", "Probe ID")
    missing_threshold = float(step_cfg.get("missing_threshold", 0.05))
    missing_code = _stringify_missing_code(step_cfg.get("missing_code", 9))

    with step_logger("Filter individuals by missingness"):
        df = pd.read_csv(input_path)
        if id_column not in df.columns:
            raise ConfigError(f"Column '{id_column}' not found in {input_path}")

        genotype_matrix = df.drop(columns=[id_column])
        genotype_str = genotype_matrix.astype(str)
        missing_fraction = genotype_str.eq(missing_code).mean(axis=0)
        keep_columns = missing_fraction[missing_fraction <= missing_threshold].index.tolist()

        filtered = pd.concat([df[[id_column]], genotype_matrix[keep_columns]], axis=1)
        ensure_parent(output_path)
        filtered.to_csv(output_path, index=False)

        if removed_log_path is not None:
            removed = missing_fraction[missing_fraction > missing_threshold]
            if not removed.empty:
                ensure_parent(removed_log_path)
                removed_df = removed.reset_index()
                removed_df.columns = ["Individual", "missing_fraction"]
                removed_df.to_csv(removed_log_path, index=False)

        logger.info(
            "Filtered individuals: kept %d of %d (missing = %.2f%%)",
            len(keep_columns),
            genotype_matrix.shape[1],
            missing_threshold * 100,
        )


def transpose_genotypes(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step04_transpose_matrix", "transpose_genotypes")
    if step_cfg is None:
        logger.info("Skipping step04_transpose_matrix (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step03_filter_individuals", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step04_transpose_matrix.output_csv must be set")

    id_column = step_cfg.get("snp_id_column", "Probe ID")
    taxa_column = step_cfg.get("taxa_column", "Taxa")

    with step_logger("Transpose genotype matrix"):
        df = pd.read_csv(input_path)
        if id_column not in df.columns:
            raise ConfigError(f"Column '{id_column}' not found in {input_path}")
        df = df.set_index(id_column)
        transposed = df.transpose().reset_index().rename(columns={"index": taxa_column})
        ensure_parent(output_path)
        transposed.to_csv(output_path, index=False)
        logger.info("Transposed matrix saved to %s", output_path)


def calculate_ibs(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step05_calculate_ibs", "ibs_similarity")
    if step_cfg is None:
        logger.info("Skipping step05_calculate_ibs (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step04_transpose_matrix", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step05_calculate_ibs.output_csv must be set")

    label_column = step_cfg.get("label_column", "Taxa")
    missing_code = step_cfg.get("missing_code", 9)

    with step_logger("Compute IBS similarity matrix"):
        df = pd.read_csv(input_path)
        if label_column not in df.columns:
            raise ConfigError(f"Column '{label_column}' not found in {input_path}")

        labels = df[label_column].astype(str).tolist()
        geno = df.drop(columns=[label_column])
        numeric = _to_numeric_dataframe(geno, missing_code)
        geno_array = numeric.to_numpy(dtype=float)

        n_samples = geno_array.shape[0]
        similarity = np.zeros((n_samples, n_samples), dtype=float)

        for i in range(n_samples):
            g1 = geno_array[i]
            for j in range(i, n_samples):
                g2 = geno_array[j]
                valid = ~np.isnan(g1) & ~np.isnan(g2)
                scored = np.count_nonzero(valid)
                if scored == 0:
                    sim = 0.0
                else:
                    g1v = g1[valid]
                    g2v = g2[valid]
                    matches = (g1v == g2v)
                    hetero = (
                        ((g1v == 1) & np.isin(g2v, (0, 2)))
                        | ((g2v == 1) & np.isin(g1v, (0, 2)))
                    )
                    common = matches.sum() + 0.5 * hetero.sum()
                    sim = common / scored
                similarity[i, j] = sim
                similarity[j, i] = sim

        matrix = pd.DataFrame(similarity, columns=labels)
        matrix.insert(0, label_column, labels)
        ensure_parent(output_path)
        matrix.to_csv(output_path, index=False)
        logger.info("IBS matrix saved to %s", output_path)


def convert_to_dissimilarity(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step06_convert_to_dissimilarity", "ibs_dissimilarity")
    if step_cfg is None:
        logger.info("Skipping step06_convert_to_dissimilarity (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step05_calculate_ibs", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step06_convert_to_dissimilarity.output_csv must be set")

    with step_logger("Convert IBS similarity to dissimilarity"):
        df = pd.read_csv(input_path)
        value_columns = df.columns[1:]
        df[value_columns] = 1 - df[value_columns]
        ensure_parent(output_path)
        df.to_csv(output_path, index=False)
        logger.info("Dissimilarity matrix saved to %s", output_path)


def extract_labels(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step07_extract_labels", "extract_labels")
    if step_cfg is None:
        logger.info("Skipping step07_extract_labels (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step05_calculate_ibs", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_txt", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step07_extract_labels.output_txt must be set")

    with step_logger("Extract genotype labels"):
        df = pd.read_csv(input_path, nrows=0)
        labels = df.columns[1:]
        ensure_parent(output_path)
        output_path.write_text("\n".join(labels) + "\n", encoding="utf-8")
        logger.info("Wrote %d labels to %s", len(labels), output_path)


def transpose_genotypes_pruned(config: PipelineConfig) -> None:
    """Same as transpose_genotypes but uses a pruned-specific config key.

    This allows running a separate flow that writes into alternative output folders
    (for example `outputs_pruned/`) without touching the original outputs.
    """
    step_cfg, step_key = _load_step_config(config, "step04_transpose_pruned", "transpose_genotypes_pruned")
    if step_cfg is None:
        logger.info("Skipping step04_transpose_pruned (not configured)")
        return

    # Reuse the same implementation but resolved from the pruned step config
    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step03_filter_individuals", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step04_transpose_pruned.output_csv must be set")

    id_column = step_cfg.get("snp_id_column", "Probe ID")
    taxa_column = step_cfg.get("taxa_column", "Taxa")

    with step_logger("Transpose genotype matrix (pruned)"):
        df = pd.read_csv(input_path)
        if id_column not in df.columns:
            raise ConfigError(f"Column '{id_column}' not found in {input_path}")
        df = df.set_index(id_column)
        transposed = df.transpose().reset_index().rename(columns={"index": taxa_column})
        ensure_parent(output_path)
        transposed.to_csv(output_path, index=False)
        logger.info("Transposed pruned matrix saved to %s", output_path)


def calculate_ibs_pruned(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step05_calculate_ibs_pruned", "ibs_similarity_pruned")
    if step_cfg is None:
        logger.info("Skipping step05_calculate_ibs_pruned (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step04_transpose_pruned", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step05_calculate_ibs_pruned.output_csv must be set")

    label_column = step_cfg.get("label_column", "Taxa")
    missing_code = step_cfg.get("missing_code", 9)

    with step_logger("Compute IBS similarity matrix (pruned)"):
        df = pd.read_csv(input_path)
        if label_column not in df.columns:
            raise ConfigError(f"Column '{label_column}' not found in {input_path}")

        labels = df[label_column].astype(str).tolist()
        geno = df.drop(columns=[label_column])
        numeric = _to_numeric_dataframe(geno, missing_code)
        geno_array = numeric.to_numpy(dtype=float)

        n_samples = geno_array.shape[0]
        similarity = np.zeros((n_samples, n_samples), dtype=float)

        for i in range(n_samples):
            g1 = geno_array[i]
            for j in range(i, n_samples):
                g2 = geno_array[j]
                valid = ~np.isnan(g1) & ~np.isnan(g2)
                scored = np.count_nonzero(valid)
                if scored == 0:
                    sim = 0.0
                else:
                    g1v = g1[valid]
                    g2v = g2[valid]
                    matches = (g1v == g2v)
                    hetero = (
                        ((g1v == 1) & np.isin(g2v, (0, 2)))
                        | ((g2v == 1) & np.isin(g1v, (0, 2)))
                    )
                    common = matches.sum() + 0.5 * hetero.sum()
                    sim = common / scored
                similarity[i, j] = sim
                similarity[j, i] = sim

        matrix = pd.DataFrame(similarity, columns=labels)
        matrix.insert(0, label_column, labels)
        ensure_parent(output_path)
        matrix.to_csv(output_path, index=False)
        logger.info("IBS matrix (pruned) saved to %s", output_path)


def convert_to_dissimilarity_pruned(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step06_convert_to_dissimilarity_pruned", "ibs_dissimilarity_pruned")
    if step_cfg is None:
        logger.info("Skipping step06_convert_to_dissimilarity_pruned (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step05_calculate_ibs_pruned", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step06_convert_to_dissimilarity_pruned.output_csv must be set")

    with step_logger("Convert IBS similarity to dissimilarity (pruned)"):
        df = pd.read_csv(input_path)
        value_columns = df.columns[1:]
        df[value_columns] = 1 - df[value_columns]
        ensure_parent(output_path)
        df.to_csv(output_path, index=False)
        logger.info("Dissimilarity matrix (pruned) saved to %s", output_path)


def extract_labels_pruned(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step07_extract_labels_pruned", "extract_labels_pruned")
    if step_cfg is None:
        logger.info("Skipping step07_extract_labels_pruned (not configured)")
        return

    input_path = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step05_calculate_ibs_pruned", "output_csv"),
    )
    output_path = config.path("steps", step_key, "output_txt", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step07_extract_labels_pruned.output_txt must be set")

    with step_logger("Extract genotype labels (pruned)"):
        df = pd.read_csv(input_path, nrows=0)
        labels = df.columns[1:]
        ensure_parent(output_path)
        output_path.write_text("\n".join(labels) + "\n", encoding="utf-8")
        logger.info("Wrote %d pruned labels to %s", len(labels), output_path)


def create_plink_map(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step08_create_plink_map", "create_plink_map")
    if step_cfg is None:
        logger.info("Skipping step08_create_plink_map (not configured)")
        return

    snp_list_path = _resolve_input_path(
        config,
        step_cfg,
        "snp_list_txt",
        fallback_keys=("steps", "step02_filter_mapped_snps", "snp_list_txt"),
    )
    annotation_path = _resolve_input_path(
        config,
        step_cfg,
        "annotation_map",
        fallback_keys=("paths", "probe_annotation_map"),
    )
    output_path = config.path("steps", step_key, "output_map", create_parent=True)
    if output_path is None:
        raise ConfigError("steps.step08_create_plink_map.output_map must be set")

    snp_id_column = step_cfg.get("snp_id_column", "35K SNPId")
    chromosome_column = step_cfg.get("chromosome_column", "IWGSC_v1_Chromosome")
    position_column = step_cfg.get("position_column", "IWGSC_v1_Position")

    with step_logger("Create PLINK map file"):
        snp_list = pd.read_csv(snp_list_path, header=None, names=["SNP"], dtype=str)
        metadata = pd.read_csv(annotation_path, sep="\t", dtype=str)
        metadata = metadata.dropna(subset=[chromosome_column, position_column])
        metadata[position_column] = pd.to_numeric(metadata[position_column], errors="coerce")
        metadata = metadata.dropna(subset=[position_column])
        metadata["Chrom"] = metadata[chromosome_column].str.replace("chr", "", regex=False)

        merged = snp_list.merge(
            metadata[[snp_id_column, "Chrom", position_column]],
            left_on="SNP",
            right_on=snp_id_column,
            how="left",
        )
        missing = merged[merged["Chrom"].isna()]
        if not missing.empty:
            logger.warning("%d SNPs missing chromosome assignments", len(missing))

        plink_map = merged[["Chrom", "SNP", position_column]].copy()
        plink_map.insert(2, "cM", 0)
        plink_map.columns = [0, 1, 2, 3]
        plink_map = plink_map.dropna()
        plink_map[0] = plink_map[0].astype(str)
        plink_map[3] = plink_map[3].astype(int)

        ensure_parent(output_path)
        plink_map.to_csv(output_path, sep="\t", header=False, index=False)
        logger.info("PLINK map saved to %s", output_path)


def transpose_to_tped(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step09_csv_to_tped", "transpose_to_tped")
    if step_cfg is None:
        logger.info("Skipping step09_csv_to_tped (not configured)")
        return

    geno_path = _resolve_input_path(
        config,
        step_cfg,
        "genotype_csv",
        fallback_keys=("steps", "step04_transpose_matrix", "output_csv"),
    )
    map_path = _resolve_input_path(
        config,
        step_cfg,
        "map_file",
        fallback_keys=("steps", "step08_create_plink_map", "output_map"),
    )
    output_tped = config.path("steps", step_key, "output_tped", create_parent=True)
    if output_tped is None:
        raise ConfigError("steps.step09_csv_to_tped.output_tped must be set")
    output_tfam = config.path("steps", step_key, "output_tfam", create_parent=True)
    if output_tfam is None:
        output_tfam = output_tped.with_suffix(".tfam")
        ensure_parent(output_tfam)

    with step_logger("Convert transposed genotype CSV to TPED"):
        geno_df = pd.read_csv(geno_path, index_col=0)
        map_df = pd.read_csv(map_path, sep="\t", header=None, names=["CHR", "SNP", "CM", "BP"])
        missing_snps = pd.Index(map_df["SNP"]).difference(geno_df.columns)
        if not missing_snps.empty:
            logger.warning("%d SNPs from map not found in genotype matrix", len(missing_snps))
        common_snps = [s for s in map_df["SNP"] if s in geno_df.columns]
        geno_df = geno_df[common_snps]

        allele_map = {
            "0": "A A",
            "1": "A B",
            "2": "B B",
            "9": "0 0",
        }

        ensure_parent(output_tped)
        with output_tped.open("w", encoding="utf-8") as handle:
            for row in map_df.itertuples(index=False):
                snp = row.SNP
                if snp not in geno_df.columns:
                    continue
                genotypes = geno_df[snp].astype(str).map(allele_map).fillna("0 0")
                line = f"{row.CHR}\t{snp}\t{row.CM}\t{row.BP}\t" + "\t".join(genotypes)
                handle.write(line + "\n")
        with output_tfam.open("w", encoding="utf-8") as tfam_handle:
            for sample in geno_df.index:
                tfam_handle.write(f"{sample}\t{sample}\t0\t0\t0\t-9\n")
        logger.info("TPED saved to %s", output_tped)
        logger.info("TFAM saved to %s", output_tfam)


def extract_pruned_subset(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step12_extract_pruned_subset", "extract_pruned_subset")
    if step_cfg is None:
        logger.info("Skipping step12_extract_pruned_subset (not configured)")
        return

    prune_in_path = _resolve_input_path(
        config,
        step_cfg,
        "prune_in",
        fallback_keys=("paths", "prune_in"),
    )
    input_csv = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
        fallback_keys=("steps", "step03_filter_individuals", "output_csv"),
    )
    output_csv = config.path("steps", step_key, "output_csv", create_parent=True)
    if output_csv is None:
        raise ConfigError("steps.step12_extract_pruned_subset.output_csv must be set")

    with step_logger("Extract SNP subset from PLINK prune.in"):
        with prune_in_path.open("r", encoding="utf-8") as handle:
            selected = {line.strip() for line in handle if line.strip()}
        with input_csv.open("r", encoding="utf-8") as fin, output_csv.open("w", encoding="utf-8") as fout:
            for idx, line in enumerate(fin):
                if idx == 0:
                    fout.write(line)
                    continue
                snp_id = line.split(",", 1)[0]
                if snp_id in selected:
                    fout.write(line)
        logger.info("Pruned subset saved to %s (%d SNPs)", output_csv, len(selected))


def convert_to_structure(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step14_convert_to_structure", "convert_to_structure")
    if step_cfg is None:
        logger.info("Skipping step14_convert_to_structure (not configured)")
        return

    input_csv = _resolve_input_path(
        config,
        step_cfg,
        "input_csv",
    fallback_keys=("steps", "step12_extract_pruned_subset", "output_csv"),
    )
    output_txt = config.path("steps", step_key, "output_txt", create_parent=True)
    if output_txt is None:
        raise ConfigError("steps.step14_convert_to_structure.output_txt must be set")

    with step_logger("Convert 012 genotypes to STRUCTURE format"):
        df = pd.read_csv(input_csv)
        if df.empty:
            raise ConfigError(f"{input_csv} has no rows to convert")

        header = df.columns[0].strip().lower()
        if header in {"probe id", "probe_id", "snp", "snp_id", "marker", "marker_id"}:
            df = df.set_index(df.columns[0]).transpose().reset_index().rename(columns={"index": "Taxa"})
        else:
            df = df.rename(columns={df.columns[0]: "Taxa"})

        encoding = {
            "0": ("1", "1"),
            "1": ("1", "2"),
            "2": ("2", "2"),
        }
        missing = ("-9", "-9")

        with output_txt.open("w", encoding="utf-8") as fout:
            for _, row in df.iterrows():
                label = str(row.iloc[0])
                alleles: List[str] = [label]
                for value in row.iloc[1:]:
                    if pd.isna(value):
                        alleles.extend(missing)
                        continue
                    code = str(value).strip()
                    if code.endswith(".0"):
                        code = code[:-2]
                    if code in encoding:
                        alleles.extend(encoding[code])
                    else:
                        alleles.extend(missing)
                fout.write("\t".join(alleles) + "\n")
        logger.info("STRUCTURE file saved to %s", output_txt)


def _resolve_rscript(config: PipelineConfig | None) -> str:
    candidate = None
    if config is not None:
        configured = config.get("tools", "rscript_executable")
        if configured:
            resolved = config.resolve_path(configured)
            candidate = str(resolved) if resolved else configured
    if not candidate:
        env_override = os.environ.get("RSCRIPT_PATH") or os.environ.get("RSCRIPT_EXECUTABLE")
        if env_override:
            candidate = env_override
    if not candidate:
        found = shutil.which("Rscript")
        if not found:
            raise ConfigError(
                "Rscript executable not found. Add it to PATH, set RSCRIPT_PATH/RSCRIPT_EXECUTABLE, "
                "or configure tools.rscript_executable."
            )
        candidate = found
    return candidate


def _run_r_script(script: Path, args: List[str], config: PipelineConfig | None = None) -> None:
    rscript = _resolve_rscript(config)
    cmd = [rscript, str(script)] + [str(arg) for arg in args]
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:  # pragma: no cover - environment dependent
        raise ConfigError("Rscript executable not found. Ensure R is installed and on PATH.") from exc
    except subprocess.CalledProcessError as exc:  # pragma: no cover - runtime error handling
        raise ConfigError(f"R script failed: {' '.join(cmd)}") from exc


def run_pcoa(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step13_pcoa", "pcoa")
    if step_cfg is None:
        logger.info("Skipping step13_pcoa (not configured)")
        return

    script = _resolve_input_path(config, step_cfg, "script")
    dissimilarity = _resolve_input_path(
        config,
        step_cfg,
        "dissimilarity_csv",
        fallback_keys=("steps", "step06_convert_to_dissimilarity", "output_csv"),
    )
    labels = _resolve_input_path(config, step_cfg, "labels_txt", fallback_keys=("steps", "step07_extract_labels", "output_txt"))
    coordinates = config.path("steps", step_key, "coordinates_csv", create_parent=True)
    plot_png = config.path("steps", step_key, "plot_png", create_parent=True)
    if coordinates is None or plot_png is None:
        raise ConfigError("pcoa coordinates_csv and plot_png must be set")

    with step_logger("Run PCoA via R"):
        _run_r_script(script, [dissimilarity, labels, coordinates, plot_png], config)
        logger.info("PCoA outputs saved to %s and %s", coordinates, plot_png)


def fix_pruned_ped(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step15_fix_pruned_ped")
    if step_cfg is None:
        logger.info("Skipping step15_fix_pruned_ped (not configured)")
        return

    input_ped = _resolve_input_path(config, step_cfg, "input_ped")
    output_ped = config.path("steps", step_key, "output_ped", create_parent=True)
    if output_ped is None:
        raise ConfigError("steps.step15_fix_pruned_ped.output_ped must be set")

    with step_logger("Normalise TASSEL PED formatting"):
        total = 0
        input_path = Path(input_ped)
        output_path = Path(output_ped)
        with input_path.open("r", encoding="utf-8") as src, output_path.open("w", encoding="utf-8") as dst:
            for line_no, raw in enumerate(src, 1):
                stripped = raw.strip()
                if not stripped:
                    continue
                tokens = stripped.split()
                if len(tokens) < 6:
                    raise ConfigError(f"Line {line_no} in {input_path} is not a valid PED record (expected ≥6 columns).")
                try:
                    sex_val = int(tokens[4])
                except ValueError as exc:
                    raise ConfigError(f"Line {line_no}: could not parse sex value '{tokens[4]}'.") from exc
                tokens[4] = str(sex_val + 1)
                for idx in range(6, len(tokens)):
                    allele = tokens[idx]
                    if allele.upper() == "B":
                        tokens[idx] = "C"
                    elif allele.upper() == "C":
                        tokens[idx] = "C"
                dst.write(" ".join(tokens) + "\n")
                total += 1
        logger.info("Wrote %s corrected PED rows to %s", total, output_path)


def split_traits(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step16_trait_split", "trait_split")
    if step_cfg is None:
        logger.info("Skipping step16_trait_split (not configured)")
        return

    script = _resolve_input_path(config, step_cfg, "script")
    mlm_stats = _resolve_input_path(config, step_cfg, "mlm_stats", fallback_keys=("paths", "mlm_stats"))
    output_dir = config.path("steps", step_key, "output_dir", create_parent=True)
    if output_dir is None:
        raise ConfigError("trait_split.output_dir must be set")

    with step_logger("Split MLM statistics by trait"):
        _run_r_script(script, [mlm_stats, output_dir], config)
        logger.info("Trait CSVs written to %s", output_dir)


def adjust_traits(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step17_trait_adjust", "trait_adjust")
    if step_cfg is None:
        logger.info("Skipping step17_trait_adjust (not configured)")
        return

    script = _resolve_input_path(config, step_cfg, "script")
    input_dir = config.path("steps", step_key, "input_dir", create_parent=False)
    if input_dir is None:
        raise ConfigError("trait_adjust.input_dir must be set")
    pattern = step_cfg.get("pattern", "^mlm_stats_.*\\.csv$")
    output_dir = config.path("steps", step_key, "output_dir", create_parent=True)
    summary_file = config.path("steps", step_key, "summary_file", create_parent=True)
    n_tests = str(step_cfg.get("n_tests", 3048))
    top_n = str(step_cfg.get("top_n", 10))

    with step_logger("Adjust trait p-values"):
        _run_r_script(script, [input_dir, pattern, output_dir, summary_file, n_tests, top_n], config)
        logger.info("Adjusted trait outputs written to %s", output_dir)


def ld_decay_plot(config: PipelineConfig) -> None:
    step_cfg, step_key = _load_step_config(config, "step10_ld_decay_plot", "ld_decay_plot")
    if step_cfg is None:
        logger.info("Skipping step10_ld_decay_plot (not configured)")
        return

    script = _resolve_input_path(config, step_cfg, "script")
    tped_path = _resolve_input_path(
        config,
        step_cfg,
        "tped",
        fallback_keys=("steps", "step09_csv_to_tped", "output_tped"),
    )
    tfam_path = _resolve_input_path(
        config,
        step_cfg,
        "tfam",
        fallback_keys=("steps", "step09_csv_to_tped", "output_tfam"),
    )
    if not tped_path.exists() or not tfam_path.exists():
        raise ConfigError("TPED/TFAM files are required for LD decay plotting (run step09_csv_to_tped first).")

    ld_summary = config.path("steps", step_key, "ld_summary", create_parent=True)
    if ld_summary is None:
        ld_summary = config.path("paths", "ld_summary", create_parent=True)
    if ld_summary is None:
        raise ConfigError("ld_decay_plot.ld_summary must be set in the configuration.")
    ensure_parent(ld_summary)

    output_png = config.path("steps", step_key, "output_png", create_parent=True)

    bins_tsv = config.path("steps", step_key, "bins_tsv", create_parent=True)
    if bins_tsv is None:
        bins_tsv = ld_summary.parent / "ld_bins.tsv"
        ensure_parent(bins_tsv)

    smooth_tsv = config.path("steps", step_key, "smooth_tsv", create_parent=True)
    if smooth_tsv is None:
        smooth_tsv = ld_summary.parent / "ld_smooth.tsv"
        ensure_parent(smooth_tsv)

    metrics_tsv = config.path("steps", step_key, "metrics_tsv", create_parent=True)
    if metrics_tsv is None:
        metrics_tsv = ld_summary.parent / "ld_decay_metrics.tsv"
        ensure_parent(metrics_tsv)

    threshold_txt = config.path("steps", step_key, "threshold_txt", create_parent=True)
    if threshold_txt is None:
        threshold_txt = ld_summary.parent / "ld_decay_threshold.txt"

    bin_size_bp = int(step_cfg.get("bin_size_bp", 50_000))
    loess_span = float(step_cfg.get("loess_span", 0.2))

    output_prefix = config.path("steps", step_key, "output_prefix", create_parent=True)
    if output_prefix is None:
        output_prefix = ld_summary.with_suffix("")
    ensure_parent(output_prefix)

    decay_cfg = config.get("tools", "plink_ld_decay", default={})
    ld_window_kb = int(step_cfg.get("ld_window_kb", decay_cfg.get("ld_window_kb", 600000)))
    ld_window = int(step_cfg.get("ld_window", decay_cfg.get("ld_window", 999999)))
    ld_window_r2 = float(step_cfg.get("ld_window_r2", decay_cfg.get("ld_window_r2", 0)))

    prune_cfg = config.get("tools", "plink_ld_prune", default={})
    r2_threshold = float(step_cfg.get("r2_threshold", prune_cfg.get("r2", 0.2)))

    explicit_plink = step_cfg.get("plink_executable") or config.get("tools", "plink_executable")
    if explicit_plink:
        plink_exec = config.resolve_path(explicit_plink) or Path(explicit_plink)
    else:
        env_plink = os.environ.get("PLINK_PATH") or os.environ.get("PLINK_EXECUTABLE")
        if env_plink:
            plink_exec = Path(env_plink)
        else:
            found = shutil.which("plink")
            if not found:
                raise ConfigError(
                    "PLINK executable not found. Add it to PATH, set PLINK_PATH, or configure tools.plink_executable."
                )
            plink_exec = Path(found)

    # Ensure Rscript is available (required for plotting)
    explicit_rscript = config.get("tools", "rscript_executable") or os.environ.get("RSCRIPT_PATH")
    if explicit_rscript:
        rscript_exec = Path(explicit_rscript)
    else:
        found_r = shutil.which("Rscript")
        if not found_r:
            raise ConfigError("Rscript not found on PATH. Install R and ensure Rscript is available or set tools.rscript_executable.")
        rscript_exec = Path(found_r)

    # Validate PLINK .map positions (4th column) are numeric before running PLINK
    # Try to locate a .map file from step08_create_plink_map or configured path
    map_path = config.path("steps", "step08_create_plink_map", "output_map")
    if map_path is not None and map_path.exists():
        try:
            with map_path.open("r", encoding="utf-8") as mf:
                bad_lines = 0
                total = 0
                for ln in mf:
                    total += 1
                    parts = ln.strip().split()
                    if len(parts) < 4:
                        bad_lines += 1
                        continue
                    try:
                        int(parts[3])
                    except ValueError:
                        bad_lines += 1
                if total > 0 and bad_lines == total:
                    raise ConfigError(f"PLINK .map file {map_path} appears to have non-numeric BP positions in column 4; check your map file.")
        except OSError:
            logger.warning("Could not open .map file %s to validate positions; continuing and letting PLINK report errors if any.", map_path)

    base_input = tped_path.with_suffix("")
    ld_prefix = Path(output_prefix)
    ld_file = ld_prefix.with_suffix(".ld")

    with step_logger("Compute LD decay with PLINK"):
        cmd = [
            str(plink_exec),
            "--tfile",
            str(base_input),
            "--allow-extra-chr",
            "--r2",
            "--ld-window", str(ld_window),
            "--ld-window-kb", str(ld_window_kb),
            "--ld-window-r2", str(ld_window_r2),
            "--out",
            str(ld_prefix),
        ]
        logger.info("Running: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise ConfigError(f"PLINK LD decay command failed: {exc.stderr or exc.stdout}") from exc

        if not ld_file.exists():
            raise ConfigError(f"Expected PLINK LD output not found: {ld_file}")

        records: list[tuple[int, float]] = []
        with ld_file.open("r", encoding="utf-8") as handle:
            header = handle.readline().strip().split()
            try:
                bp_a_idx = header.index("BP_A")
                bp_b_idx = header.index("BP_B")
                r2_idx = header.index("R2")
            except ValueError as exc:
                raise ConfigError("Unexpected PLINK LD output format; expected columns BP_A, BP_B, R2.") from exc
            for line in handle:
                if not line.strip():
                    continue
                parts = line.split()
                try:
                    bp_a = int(parts[bp_a_idx])
                    bp_b = int(parts[bp_b_idx])
                    r2_val = float(parts[r2_idx])
                except (ValueError, IndexError):
                    continue
                records.append((abs(bp_b - bp_a), r2_val))

        if not records:
            raise ConfigError("No LD records were produced; adjust ld_window parameters or verify the dataset.")

        records.sort(key=lambda item: item[0])
        with ld_summary.open("w", encoding="utf-8") as summary_handle:
            for dist, r2_val in records:
                summary_handle.write(f"{dist}\t{r2_val:.6f}\n")
        _run_r_script(
            script,
            [
                ld_summary,
                bins_tsv,
                smooth_tsv,
                metrics_tsv,
                output_png,
                str(int(bin_size_bp)),
                f"{loess_span}",
            ],
            config,
        )
        logger.info("LD decay analysis outputs saved to %s and %s", output_png, metrics_tsv)

    if not metrics_tsv.exists():
        raise ConfigError(f"Expected LD metrics file not found: {metrics_tsv}")

    try:
        metrics_df = pd.read_csv(metrics_tsv, sep="\t", dtype={"metric": str, "value": object})
    except Exception as exc:  # pragma: no cover - defensive
        raise ConfigError(f"Failed to read LD metrics TSV {metrics_tsv}: {exc}") from exc

    metric_map = {}
    for _, row in metrics_df.iterrows():
        metric = str(row.get("metric", "")).strip()
        value = row.get("value")
        if metric:
            metric_map[metric] = value

    def _metric_to_float(name: str) -> float | None:
        value = metric_map.get(name)
        if pd.isna(value):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    if metric_map.get("warnings"):
        warnings_text = str(metric_map["warnings"]).strip()
        if warnings_text:
            logger.warning("LD decay warnings: %s", warnings_text)

    d0_bin = _metric_to_float("d0.2_bin_bp")
    d0_loess = _metric_to_float("d0.2_loess_bp")

    ensure_parent(threshold_txt)
    with threshold_txt.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"r2_threshold\t{r2_threshold}\n")
        if d0_bin is not None:
            handle.write(f"d0.2_bin_bp\t{int(round(d0_bin))}\n")
            handle.write(f"d0.2_bin_kb\t{d0_bin / 1000.0:.2f}\n")
        else:
            handle.write("d0.2_bin_bp\tNA\n")
            handle.write("d0.2_bin_kb\tNA\n")
        if d0_loess is not None:
            handle.write(f"d0.2_loess_bp\t{int(round(d0_loess))}\n")
            handle.write(f"d0.2_loess_kb\t{d0_loess / 1000.0:.2f}\n")
        else:
            handle.write("d0.2_loess_bp\tNA\n")
            handle.write("d0.2_loess_kb\tNA\n")
    logger.info("LD decay metrics written to %s", threshold_txt)




STEP_FUNCTIONS = {
    "step01_filter_variants": filter_variants,
    "step02_filter_mapped_snps": filter_mapped_snps,
    "step03_filter_individuals": filter_individuals,
    "step04_transpose_matrix": transpose_genotypes,
    "step05_calculate_ibs": calculate_ibs,
    "step06_convert_to_dissimilarity": convert_to_dissimilarity,
    "step07_extract_labels": extract_labels,
    "step08_create_plink_map": create_plink_map,
    "step09_csv_to_tped": transpose_to_tped,
    "step10_ld_decay_plot": ld_decay_plot,
    "step11_plink_pruning": plink_pruning,
    "step12_extract_pruned_subset": extract_pruned_subset,
    "step13_pcoa": run_pcoa,
    "step14_convert_to_structure": convert_to_structure,
    "step15_fix_pruned_ped": fix_pruned_ped,
    "step16_trait_split": split_traits,
    "step17_trait_adjust": adjust_traits,
    "step04_transpose_pruned": transpose_genotypes_pruned,
    "step05_calculate_ibs_pruned": calculate_ibs_pruned,
    "step06_convert_to_dissimilarity_pruned": convert_to_dissimilarity_pruned,
    "step07_extract_labels_pruned": extract_labels_pruned,
}

# Backward compatibility with legacy step names
STEP_FUNCTIONS.update({
    "filter_variants": filter_variants,
    "filter_mapped_snps": filter_mapped_snps,
    "filter_individuals": filter_individuals,
    "transpose_genotypes": transpose_genotypes,
    "ibs_similarity": calculate_ibs,
    "ibs_dissimilarity": convert_to_dissimilarity,
    "extract_labels": extract_labels,
    "create_plink_map": create_plink_map,
    "transpose_to_tped": transpose_to_tped,
    "plink_pruning": plink_pruning,
    "extract_pruned_subset": extract_pruned_subset,
    "convert_to_structure": convert_to_structure,
    "pcoa": run_pcoa,
    "fix_pruned_ped": fix_pruned_ped,
    "trait_split": split_traits,
    "trait_adjust": adjust_traits,
    "ld_decay_plot": ld_decay_plot,
})




