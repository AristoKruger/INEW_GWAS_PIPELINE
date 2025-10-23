from __future__ import annotations

import argparse
import sys

from .config import ConfigError, load_config
from .pipeline import GWASPipeline
from .utils import get_logger


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="INEW GWAS pipeline runner")
    parser.add_argument(
        "--config",
        default="config/pipeline.yaml",
        help="Path to the pipeline YAML configuration file",
    )
    parser.add_argument(
        "--list-steps",
        action="store_true",
        help="List all available steps and exit",
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        help="Subset of steps to run (default: run all configured steps)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the steps that would run without executing them",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logger = get_logger()

    try:
        config = load_config(args.config)
    except ConfigError as exc:
        logger.error(str(exc))
        return 1

    pipeline = GWASPipeline(config)

    if args.list_steps:
        logger.info("Available steps: %s", ", ".join(pipeline.available_steps()))
        logger.info("Configured steps: %s", ", ".join(pipeline.configured_steps()))
        return 0

    try:
        pipeline.run(selected_steps=args.steps, dry_run=args.dry_run)
    except ConfigError as exc:
        logger.error(str(exc))
        return 1

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
