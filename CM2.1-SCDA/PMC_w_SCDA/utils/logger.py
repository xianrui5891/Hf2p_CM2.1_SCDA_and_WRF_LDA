from __future__ import annotations

import logging
import os
import sys
from datetime import datetime
from pathlib import Path


_GLOBAL_LOGGER: logging.LoggerAdapter | None = None


def _rank() -> int:
    for name in ("RANK", "LOCAL_RANK", "OMPI_COMM_WORLD_RANK", "PMI_RANK"):
        if name in os.environ:
            try:
                return int(os.environ[name])
            except ValueError:
                pass
    return 0


class _RankFilter(logging.Filter):
    def __init__(self, rank: int) -> None:
        super().__init__()
        self.rank = rank

    def filter(self, record: logging.LogRecord) -> bool:
        record.rank = self.rank
        return True


def configure_logger(
    run_name: str = "cm2_ae",
    output_dir: str | Path | None = None,
    console_level: str | int = "INFO",
    file_level: str | int | None = None,
) -> logging.LoggerAdapter:
    global _GLOBAL_LOGGER
    rank = _rank()
    logger = logging.getLogger(f"cm2_lda_rank{rank}")
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)
    logger.propagate = False
    logger.setLevel(logging.DEBUG)
    logger.addFilter(_RankFilter(rank))

    if isinstance(console_level, str):
        console_level = getattr(logging, console_level.upper(), logging.INFO)
    if isinstance(file_level, str):
        file_level = getattr(logging, file_level.upper(), logging.DEBUG)

    fmt = logging.Formatter("%(asctime)s [rank %(rank)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(console_level if rank == 0 else logging.WARNING)
    console.setFormatter(fmt)
    logger.addHandler(console)

    run_dir = None
    if output_dir and file_level is not None:
        run_dir = Path(output_dir) / run_name
        run_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(run_dir / f"rank{rank}.log", mode="a", encoding="utf-8")
        file_handler.setLevel(file_level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    adapter = logging.LoggerAdapter(logger, {"rank": rank})
    if run_dir is None:
        adapter.info("Logging initialized on console only; capture with nohup redirection.")
    else:
        adapter.info("Logging initialized in %s", run_dir)
    _GLOBAL_LOGGER = adapter
    return adapter


def get_logger() -> logging.LoggerAdapter:
    global _GLOBAL_LOGGER
    if _GLOBAL_LOGGER is None:
        return configure_logger(run_name=f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    return _GLOBAL_LOGGER
