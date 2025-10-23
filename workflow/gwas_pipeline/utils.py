from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path
from time import perf_counter


LOGGER_NAME = "gwas_pipeline"


def get_logger() -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%H:%M:%S")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


@contextmanager
def step_logger(name: str):
    logger = get_logger()
    logger.info("? %s", name)
    start = perf_counter()
    try:
        yield logger
    finally:
        duration = perf_counter() - start
        logger.info("? %s (%.2fs)", name, duration)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, content: str) -> None:
    ensure_parent(path)
    path.write_text(content, encoding="utf-8")
