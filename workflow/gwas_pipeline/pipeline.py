from __future__ import annotations

from datetime import datetime
import logging
from typing import Iterable, List, Optional, Sequence

from .config import ConfigError, PipelineConfig
from .progress import ProgressLog
from .steps import STEP_FUNCTIONS
from .utils import get_logger


STEP_SEQUENCE: Sequence[tuple[str, str]] = (
    ("step01_filter_variants", "filter_variants"),
    ("step02_filter_mapped_snps", "filter_mapped_snps"),
    ("step03_filter_individuals", "filter_individuals"),
    ("step04_transpose_matrix", "transpose_genotypes"),
    ("step05_calculate_ibs", "ibs_similarity"),
    ("step06_convert_to_dissimilarity", "ibs_dissimilarity"),
    ("step07_extract_labels", "extract_labels"),
    ("step08_create_plink_map", "create_plink_map"),
    ("step09_csv_to_tped", "transpose_to_tped"),
    ("step10_ld_decay_plot", "ld_decay_plot"),
    ("step11_plink_pruning", "plink_pruning"),
    ("step12_extract_pruned_subset", "extract_pruned_subset"),
    ("step13_pcoa", "pcoa"),
    ("step14_convert_to_structure", "convert_to_structure"),
    ("step15_fix_pruned_ped", "fix_pruned_ped"),
    ("step16_trait_split", "trait_split"),
    ("step17_trait_adjust", "trait_adjust"),
)

STEP_ORDER: Sequence[str] = tuple(primary for primary, _ in STEP_SEQUENCE)
STEP_ALIASES = {alias: primary for primary, alias in STEP_SEQUENCE}
STEP_ALIASES.update({primary: primary for primary, _ in STEP_SEQUENCE})


class GWASPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = get_logger()

    def configured_steps(self) -> List[str]:
        steps_cfg = self.config.get("steps", default={})
        if not isinstance(steps_cfg, dict):
            raise ConfigError("steps must be a mapping in the YAML configuration")
        configured: List[str] = []
        for primary, legacy in STEP_SEQUENCE:
            cfg = steps_cfg.get(primary)
            if not isinstance(cfg, dict) and legacy:
                cfg = steps_cfg.get(legacy)
            if isinstance(cfg, dict):
                configured.append(primary)
        return configured

    def available_steps(self) -> List[str]:
        return sorted(STEP_FUNCTIONS.keys())

    def run(self, selected_steps: Optional[Iterable[str]] = None, dry_run: bool = False) -> None:
        if selected_steps is not None:
            selected_list: List[str] = []
            for token in selected_steps:
                if token is None:
                    continue
                parts = str(token).split(",")
                for part in parts:
                    part = part.strip()
                    if part:
                        selected_list.append(part)
            normalised: List[str] = []
            unknown: List[str] = []
            for step in selected_list:
                resolved = STEP_ALIASES.get(step, step)
                if resolved in STEP_FUNCTIONS:
                    normalised.append(resolved)
                else:
                    unknown.append(step)
            if unknown:
                raise ConfigError(f"Unknown step(s): {', '.join(unknown)}")

            run_order: List[str] = []
            seen: set[str] = set()
            for step in normalised:
                if step in seen:
                    continue
                seen.add(step)
                run_order.append(step)
            if not run_order:
                self.logger.warning("No matching steps to run for selection: %s", selected_list)
                return
        else:
            run_order = self.configured_steps()

        log_handler: logging.Handler | None = None
        progress: ProgressLog | None = None
        current_step: Optional[str] = None
        try:
            log_dir = self.config.root / "outputs" / "run_history"
            log_dir.mkdir(parents=True, exist_ok=True)
            start_time = datetime.now()
            suffix = "_dryrun" if dry_run else ""
            session_name = f"run_{start_time.strftime('%Y%m%d_%H%M%S')}{suffix}"
            log_path = log_dir / f"{session_name}.log"
            progress_path = log_dir / f"{session_name}.json"
            progress = ProgressLog(
                progress_path,
                session=session_name,
                dry_run=dry_run,
                run_order=run_order,
                started_at=start_time,
            )

            log_handler = logging.FileHandler(log_path, encoding="utf-8")
            log_handler.setFormatter(
                logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%H:%M:%S")
            )
            self.logger.addHandler(log_handler)
            self.logger.info("Run started (dry_run=%s): %s", dry_run, ", ".join(run_order))

            if dry_run:
                for step in run_order:
                    progress.plan_step(step)
                    self.logger.info("[dry-run] %s", step)
                progress.finish(status="dry_run")
                self.logger.info("Dry run complete")
                return

            for step in run_order:
                current_step = step
                progress.start_step(step)
                handler = STEP_FUNCTIONS[step]
                handler(self.config)
                progress.complete_step(step)
            progress.finish(status="completed")
            self.logger.info("Run complete")
        except Exception as exc:  # pragma: no cover - defensive logging
            if progress is not None:
                if current_step is not None:
                    try:
                        progress.fail_step(current_step, message=str(exc))
                    except Exception:
                        # Best effort: avoid shadowing the original error.
                        pass
                try:
                    progress.finish(status="failed", message=str(exc))
                except Exception:
                    pass
            raise
        finally:
            if log_handler is not None:
                self.logger.removeHandler(log_handler)
                log_handler.close()


__all__ = ["GWASPipeline", "STEP_ORDER"]
