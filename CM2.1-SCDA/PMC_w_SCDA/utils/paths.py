from __future__ import annotations

import json
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = PROJECT_DIR.parent
TEMP_DIR = ROOT_DIR / "temp"
NMC_TESTER_DIR = ROOT_DIR / "nmc_tester"
PATH_CONFIG = PROJECT_DIR / "configs" / "paths.json"


def _load_path_config() -> dict[str, str]:
    if not PATH_CONFIG.exists():
        return {}
    with open(PATH_CONFIG, "r", encoding="utf-8") as f:
        return json.load(f)


def _as_path(value: str | Path, fallback: str | Path) -> Path:
    path = Path(value or fallback)
    if path.is_absolute():
        return path
    return (PROJECT_DIR / path).resolve()


_PATHS = _load_path_config()

STORAGE_ROOT = _as_path(_PATHS.get("storage_root"), "/data/cm2_lda")
RAW_DIR = _as_path(_PATHS.get("raw_dir"), STORAGE_ROOT / "raw")
DATA_DIR = _as_path(_PATHS.get("data_dir"), STORAGE_ROOT / "data")
RUNS_DIR = _as_path(_PATHS.get("runs_dir"), STORAGE_ROOT / "checkpoints")
ANALYSIS_DIR = _as_path(_PATHS.get("analysis_dir"), STORAGE_ROOT / "analysis_outputs")
NMC_OUTPUT_DIR = _as_path(_PATHS.get("nmc_output_dir"), STORAGE_ROOT / "nmc_outputs")

DEFAULT_PROCESSED_PT = _as_path(_PATHS.get("processed_data_path"), DATA_DIR)
DEFAULT_METADATA_PATH = _as_path(_PATHS.get("metadata_path"), DATA_DIR / "normalization_metadata.json")
DEFAULT_NORMALIZATION_STATS_PATH = _as_path(_PATHS.get("normalization_stats_path"), DATA_DIR / "normalization_stats.json")
DEFAULT_OBS_PATH = _as_path(_PATHS.get("observation_path"), DATA_DIR / "merged_observation.pt")
DEFAULT_AE_MODEL_PATH = _as_path(_PATHS.get("ae_model_path"), RUNS_DIR / "run_ae" / "best_model.pth")
DEFAULT_NMC_BACKGROUND_COVARIANCE_PATH = str(
    _as_path(_PATHS.get("nmc_background_covariance_path"), NMC_OUTPUT_DIR / "latent_nmc_background_covariance.npz")
)


def ensure_project_dirs() -> None:
    for path in (RAW_DIR, DATA_DIR, RUNS_DIR, ANALYSIS_DIR, NMC_OUTPUT_DIR):
        path.mkdir(parents=True, exist_ok=True)
