#!/usr/bin/env python
"""Quick health check for pipeline inputs and generated folders."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import yaml


OK_MARK = "[OK]"
MISS_MARK = "[MISSING]"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate expected inputs for the GWAS pipeline")
    parser.add_argument(
        "--config",
        default="config/pipeline.yaml",
        type=Path,
        help="Path to the pipeline configuration file",
    )
    return parser


def load_config(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"Config file not found: {path}")
    return yaml.safe_load(path.read_text()) or {}


def check_files(paths: Iterable[tuple[str, Path]]) -> bool:
    ok = True
    for label, path in paths:
        if path and path.exists():
            print(f"{OK_MARK} {label}: {path}")
        else:
            print(f"{MISS_MARK} {label}: {path if path else 'not set'}")
            ok = False
    return ok


def main() -> int:
    args = build_parser().parse_args()
    config = load_config(args.config)

    paths = config.get("paths", {})
    root = args.config.parent.parent

    expected_files = [
        ("Raw genotypes", root / paths.get("raw_genotypes", "")),
        ("Probe annotation (Excel)", root / paths.get("probe_annotation_excel", "")),
        ("Probe annotation (map)", root / paths.get("probe_annotation_map", "")),
    ]

    print("Checking required input files...")
    inputs_ok = check_files(expected_files)

    generated_files = [
        ("PLINK prune list", root / paths.get("prune_in", "")),
        ("TASSEL MLM stats", root / paths.get("mlm_stats", "")),
        ("LD summary", root / paths.get("ld_summary", "")),
    ]

    print("\nChecking generated outputs (ignore missing files before you run those steps)...")
    check_files(generated_files)

    return 0 if inputs_ok else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
