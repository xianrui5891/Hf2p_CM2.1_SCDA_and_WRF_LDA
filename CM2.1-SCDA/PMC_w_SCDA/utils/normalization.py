from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np
import torch


DEFAULT_INVALID_ABS_THRESHOLD = float("inf")


def _to_float(value: object, default: float) -> float:
    try:
        return float(value)
    except Exception:
        return float(default)


class NormalizationConfig:
    """PALM-GDA style z-score normalization, with metadata fallback."""

    def __init__(self, metadata_path: str | Path | None = None, metadata_dict: Mapping | None = None) -> None:
        base_dir: Path | None = None
        if metadata_dict is not None:
            metadata = dict(metadata_dict)
        elif metadata_path is not None:
            metadata_path = Path(metadata_path)
            base_dir = metadata_path.parent
            with open(metadata_path, "r", encoding="utf-8") as f:
                metadata = json.load(f)
        else:
            raise ValueError("metadata_path or metadata_dict is required")

        self.metadata = metadata
        self.atm_vars = list(metadata.get("atm_vars", ["TEMP", "UCOMP", "VCOMP", "PS"]))
        self.ocn_vars = list(metadata.get("ocn_vars", ["SST", "U_SURF", "V_SURF", "ETA_T"]))
        self.handle_nan = bool(metadata.get("handle_nan", True))
        self.nan_fill_value = _to_float(metadata.get("nan_fill_value", 0.0), 0.0)
        self.invalid_abs_threshold = _to_float(
            metadata.get("invalid_abs_threshold", DEFAULT_INVALID_ABS_THRESHOLD),
            DEFAULT_INVALID_ABS_THRESHOLD,
        )

        self.stats = self._load_stats(metadata, base_dir)

    def _load_stats(self, metadata: Mapping, base_dir: Path | None = None) -> dict[str, dict[str, float]]:
        raw = metadata.get("normalization_stats") or metadata.get("norm_stats")
        if raw is None:
            stats_path = metadata.get("normalization_stats_path") or metadata.get("normalization_stats_file")
            if stats_path:
                path = Path(str(stats_path))
                if not path.is_absolute() and base_dir is not None:
                    path = base_dir / path
                if path.exists():
                    with open(path, "r", encoding="utf-8") as f:
                        raw = json.load(f)
        stats: dict[str, dict[str, float]] = {}
        if isinstance(raw, Mapping):
            for name, values in raw.items():
                if isinstance(values, Mapping):
                    stats[str(name)] = {
                        "mean": _to_float(values.get("mean", 0.0), 0.0),
                        "std": max(_to_float(values.get("std", 1.0), 1.0), 1e-8),
                    }

        global_min = metadata.get("global_min", {})
        global_max = metadata.get("global_max", {})
        if isinstance(global_min, Mapping) and isinstance(global_max, Mapping):
            for key, min_value in global_min.items():
                if key in stats or key not in global_max:
                    continue
                mn = _to_float(min_value, 0.0)
                mx = _to_float(global_max[key], mn + 1.0)
                stats[str(key)] = {
                    "mean": 0.5 * (mn + mx),
                    "std": max((mx - mn) / 6.0, 1e-8),
                }
        return stats

    def get_var_stat(self, var_name: str, domain: str) -> tuple[float, float]:
        candidates = (f"{domain}::{var_name}", var_name, var_name.upper())
        for key in candidates:
            if key in self.stats:
                stat = self.stats[key]
                return stat["mean"], max(stat["std"], 1e-8)
        raise KeyError(f"Normalization stats for {domain}::{var_name} are missing.")


def normalize_standardize(
    data: np.ndarray | torch.Tensor,
    mean: float,
    std: float,
    handle_nan: bool = True,
    nan_fill_value: float = 0.0,
    invalid_abs_threshold: float = DEFAULT_INVALID_ABS_THRESHOLD,
) -> np.ndarray | torch.Tensor:
    is_torch = isinstance(data, torch.Tensor)
    if handle_nan:
        if is_torch:
            mask = torch.isnan(data) | torch.isinf(data) | (data.abs() > invalid_abs_threshold)
            clean = data.clone()
            clean[mask] = nan_fill_value
        else:
            mask = np.isnan(data) | np.isinf(data) | (np.abs(data) > invalid_abs_threshold)
            clean = np.array(data, copy=True)
            clean[mask] = nan_fill_value
    else:
        clean = data
    return (clean - mean) / max(float(std), 1e-8)


def denormalize_standardize(data: np.ndarray | torch.Tensor, mean: float, std: float) -> np.ndarray | torch.Tensor:
    return data * max(float(std), 1e-8) + mean


class VariableNormalizer:
    def __init__(self, config: NormalizationConfig) -> None:
        self.config = config

    def normalize(
        self,
        data: np.ndarray | torch.Tensor,
        var_name: str,
        domain: str = "atm",
        method: str = "standardize",
        target_range: tuple[float, float] | None = None,
    ) -> np.ndarray | torch.Tensor:
        del target_range
        if method not in {"standardize", "zscore", "palm"}:
            raise ValueError("CM2-lda uses PALM-GDA z-score normalization; use method='standardize'.")
        mean, std = self.config.get_var_stat(var_name, domain)
        return normalize_standardize(
            data,
            mean,
            std,
            handle_nan=self.config.handle_nan,
            nan_fill_value=self.config.nan_fill_value,
            invalid_abs_threshold=self.config.invalid_abs_threshold,
        )

    def denormalize(
        self,
        normalized_data: np.ndarray | torch.Tensor,
        var_name: str,
        domain: str = "atm",
        method: str = "standardize",
        source_range: tuple[float, float] | None = None,
    ) -> np.ndarray | torch.Tensor:
        del source_range
        if method not in {"standardize", "zscore", "palm"}:
            raise ValueError("CM2-lda uses PALM-GDA z-score normalization; use method='standardize'.")
        mean, std = self.config.get_var_stat(var_name, domain)
        return denormalize_standardize(normalized_data, mean, std)
