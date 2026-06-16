from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import torch

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from data_preprocess.nmc_background import (  # noqa: E402
    DEFAULT_MAX_SAMPLES,
    OnlineNMCStats,
    _autocast_context,
    _optional_int_list,
    _unpack_atm_ocn_batch,
    build_loader,
    compute_diag_approx_diagnostics,
    compute_offdiag_correlation_distribution,
    compute_sampled_diag_approx_diagnostics,
    configure_logging,
    load_model,
    plot_bdiag,
    plot_correlation_matrix,
    plot_diagonal_band_mosaic,
    plot_diagonal_block_mosaic,
    plot_diagonal_block_tiles,
    plot_local_covariance_block,
    plot_matrix,
    plot_offdiag_correlation_distribution,
    plot_paper_correlation_summary,
    resolve_device,
)
from utils.data import CM2DatasetView, CM2FrameDataset, CM2TensorDataset, FRAME_MANIFEST  # noqa: E402

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None


LOGGER = logging.getLogger("cm2_lda.climatological_background")


def _safe_name(value: str) -> str:
    name = Path(value).stem
    name = re.sub(r"^checkpoint_", "ckpt_", name)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return name or "clim_version"


def _resolve_project_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (PROJECT_DIR / path).resolve()


def _parse_datetime(value: Any, index: int, config: dict[str, Any]) -> datetime | None:
    if isinstance(value, datetime):
        return value
    if isinstance(value, str):
        text = value.strip()
        if text:
            try:
                return datetime.fromisoformat(text.replace("Z", "+00:00")).replace(tzinfo=None)
            except ValueError:
                for fmt in ("%Y-%m-%d_%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y%m%d%H", "%Y%m%d"):
                    try:
                        return datetime.strptime(text, fmt)
                    except ValueError:
                        pass
    try:
        numeric = float(value)
    except Exception:
        numeric = math.nan
    base_time = config.get("base_time")
    if base_time and math.isfinite(numeric):
        base = datetime.fromisoformat(str(base_time).replace("Z", "+00:00")).replace(tzinfo=None)
        units = str(config.get("time_value_units", "hours")).lower()
        if units in {"day", "days"}:
            return base + timedelta(days=numeric)
        if units in {"hour", "hours"}:
            return base + timedelta(hours=numeric)
        if units in {"second", "seconds"}:
            return base + timedelta(seconds=numeric)
        if units in {"index", "frame", "frames", "step", "steps"}:
            interval = float(config.get("sample_interval_hours", 1.0))
            return base + timedelta(hours=numeric * interval)
    if config.get("base_time"):
        base = datetime.fromisoformat(str(config["base_time"]).replace("Z", "+00:00")).replace(tzinfo=None)
        interval = float(config.get("sample_interval_hours", 1.0))
        return base + timedelta(hours=index * interval)
    return None


def _climatology_key(value: Any, index: int, config: dict[str, Any]) -> tuple[int, ...]:
    group = str(config.get("climatology_group", "month")).lower()
    dt = _parse_datetime(value, index, config)
    if dt is None:
        samples_per_year = int(config.get("samples_per_year", 12))
        samples_per_group = max(int(round(samples_per_year / 12)), 1)
        month = (index // samples_per_group) % 12 + 1
        hour = int(config.get("fallback_hour", 0))
    else:
        month = int(dt.month)
        hour = int(dt.hour)
    if group in {"month_hour", "month-hour", "monthhour", "diurnal", "monthly_diurnal"}:
        return month, hour
    if group in {"month", "monthly"}:
        return (month,)
    if group in {"all", "global", "none"}:
        return (0,)
    raise ValueError(f"Unsupported climatology_group={group!r}. Use month, month_hour, or all.")


def _build_dataset(config: dict[str, Any]):
    data_path = _resolve_project_path(config["data_path"])
    max_samples = config.get("max_samples")
    max_samples = None if max_samples is None else int(max_samples)
    split = str(config.get("split", "all")).lower()
    val_fraction = float(config.get("val_fraction", 0.1))
    if data_path.is_dir() and (data_path / FRAME_MANIFEST).exists():
        return CM2FrameDataset(
            data_path,
            split=split,
            val_fraction=val_fraction,
            max_samples=max_samples,
            load_mmap=bool(config.get("mmap_frames", False)),
        )
    if data_path.is_dir():
        raise FileNotFoundError(f"{FRAME_MANIFEST} not found in {data_path}")
    dataset = CM2TensorDataset(data_path, metadata_path=config.get("metadata_path"), normalize=True)
    if split != "all":
        total = len(dataset)
        n_val = max(1, int(total * val_fraction)) if total > 1 else 0
        order = list(range(total))
        indices = order[-n_val:] if split == "val" and n_val > 0 else order[:-n_val] if n_val > 0 else order
        dataset = CM2DatasetView(dataset, indices)
    if max_samples is not None:
        dataset = CM2DatasetView(dataset, range(min(max_samples, len(dataset))))
    return dataset


def _time_values(dataset, sample_count: int) -> list[Any]:
    values = list(getattr(dataset, "time_values", list(range(len(dataset)))))
    return values[:sample_count]


@torch.inference_mode()
def _climatology_pass(
    model,
    dataset,
    device: torch.device,
    batch_size: int,
    max_samples: int,
    config: dict[str, Any],
) -> tuple[dict[tuple[int, ...], np.ndarray], dict[tuple[int, ...], int], tuple[int, ...], int, list[Any]]:
    loader = build_loader(dataset, batch_size=batch_size, config=config)
    total_samples = min(len(dataset), max_samples)
    time_values = _time_values(dataset, total_samples)
    group_sum: dict[tuple[int, ...], np.ndarray] = {}
    group_count: dict[tuple[int, ...], int] = {}
    latent_shape: tuple[int, ...] | None = None
    seen = 0
    progress = tqdm(total=total_samples, desc="Pass 1/2 latent climatology", unit="sample") if tqdm is not None else None
    for batch in loader:
        if seen >= total_samples:
            break
        atm, ocn = _unpack_atm_ocn_batch(batch)
        keep = min(int(atm.shape[0]), total_samples - seen)
        if keep <= 0:
            break
        atm = torch.nan_to_num(atm[:keep].to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        ocn = torch.nan_to_num(ocn[:keep].to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        with _autocast_context(device, config):
            latent = model.encode(atm, ocn)
        if latent_shape is None:
            latent_shape = tuple(latent.shape[1:])
            LOGGER.info("Latent shape=%s latent_dim=%d", latent_shape, int(np.prod(latent_shape)))
        flat = latent.reshape(latent.shape[0], -1).float().cpu().numpy()
        for local_index, row in enumerate(flat):
            sample_index = seen + local_index
            key = _climatology_key(time_values[sample_index], sample_index, config)
            if key not in group_sum:
                group_sum[key] = np.zeros(row.shape[0], dtype=np.float64)
                group_count[key] = 0
            group_sum[key] += row.astype(np.float64, copy=False)
            group_count[key] += 1
        seen += keep
        if progress is not None:
            progress.update(keep)
    if progress is not None:
        progress.close()
    if latent_shape is None:
        raise RuntimeError("No samples were encoded for climatology statistics.")
    climatology = {
        key: (group_sum[key] / max(group_count[key], 1)).astype(np.float32)
        for key in sorted(group_sum)
    }
    return climatology, group_count, latent_shape, seen, time_values


def _make_stats(latent_shape: tuple[int, ...], plot_dims: int, config: dict[str, Any]) -> OnlineNMCStats:
    latent_dim = int(np.prod(latent_shape))
    return OnlineNMCStats(
        latent_dim=latent_dim,
        latent_shape=latent_shape,
        local_dims=plot_dims,
        diagnostic_pair_count=int(config.get("diagnostic_pair_count", 1_000_000)),
        diagnostic_seed=int(config.get("diagnostic_seed", 20260525)),
        pair_chunk_size=int(config.get("diagnostic_pair_chunk_size", 200_000)),
        diagonal_band_radius=int(config.get("diagonal_band_radius", 0)),
        diagonal_block_size=int(config.get("diagonal_block_size", config.get("plot_dims", 1))),
        diagonal_block_count=int(config.get("diagonal_block_count", 0)),
        diagonal_block_starts=_optional_int_list(config.get("diagonal_block_starts")),
        channel_spatial_diagnostics=bool(config.get("channel_spatial_diagnostics", False)),
        spatial_channel_index=int(config.get("spatial_channel_index", 0)),
        all_spatial_channel_correlations=bool(config.get("plot_all_spatial_channel_correlations", False)),
    )


def _scale_covariance_outputs(stats: dict[str, Any], alpha: float) -> dict[str, Any]:
    alpha2 = float(alpha) ** 2
    if alpha2 == 1.0:
        return stats
    for key in (
        "variance",
        "B_diag",
        "local_covariance",
        "diagonal_band_covariance",
        "diagonal_block_covariance",
        "sampled_offdiag_covariance",
    ):
        if key in stats:
            stats[key] = (np.asarray(stats[key]) * alpha2).astype(np.asarray(stats[key]).dtype, copy=False)
    return stats


@torch.inference_mode()
def encode_and_compute_climatological_bz(
    model,
    dataset,
    device: torch.device,
    batch_size: int,
    max_samples: int,
    plot_dims: int,
    config: dict[str, Any],
) -> tuple[dict[str, Any], tuple[int, ...], int, list[Any], dict[tuple[int, ...], int]]:
    start_time = time.perf_counter()
    LOGGER.info(
        "Climatological Bz start: samples=%d group=%s batch_size=%d loader_workers=%d",
        min(len(dataset), max_samples),
        config.get("climatology_group", "month"),
        batch_size,
        int(config.get("num_workers", min(4, os.cpu_count() or 4))),
    )
    climatology, group_count, latent_shape, encoded_samples, time_values = _climatology_pass(
        model,
        dataset,
        device=device,
        batch_size=batch_size,
        max_samples=max_samples,
        config=config,
    )
    LOGGER.info("Climatology groups: %s", {str(key): group_count[key] for key in sorted(group_count)})

    loader = build_loader(dataset, batch_size=batch_size, config=config)
    stats = _make_stats(latent_shape, plot_dims, config)
    seen = 0
    log_every_batches = int(config.get("log_every_n_batches", 200))
    progress = tqdm(total=encoded_samples, desc="Pass 2/2 latent anomalies", unit="sample") if tqdm is not None else None
    for batch_index, batch in enumerate(loader, start=1):
        if seen >= encoded_samples:
            break
        atm, ocn = _unpack_atm_ocn_batch(batch)
        keep = min(int(atm.shape[0]), encoded_samples - seen)
        if keep <= 0:
            break
        atm = torch.nan_to_num(atm[:keep].to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        ocn = torch.nan_to_num(ocn[:keep].to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        with _autocast_context(device, config):
            latent = model.encode(atm, ocn)
        flat = latent.reshape(latent.shape[0], -1).float().cpu().numpy()
        anomalies = np.empty_like(flat, dtype=np.float32)
        for local_index, row in enumerate(flat):
            sample_index = seen + local_index
            key = _climatology_key(time_values[sample_index], sample_index, config)
            anomalies[local_index] = row - climatology[key]
        stats.update(anomalies)
        seen += keep
        if progress is not None:
            progress.update(keep)
        if batch_index % log_every_batches == 0:
            elapsed = max(time.perf_counter() - start_time, 1e-6)
            LOGGER.info("Anomalies encoded %d/%d, rate=%.2f sample/s", seen, encoded_samples, seen / elapsed)
    if progress is not None:
        progress.close()
    final = stats.finalize()
    final["sample_count"] = int(seen)
    final["valid_pairs"] = int(final["valid_pairs"])
    final = _scale_covariance_outputs(final, float(config.get("alpha", 1.0)))
    LOGGER.info("Climatological Bz done: samples=%d anomalies=%d elapsed=%.1fs", seen, int(final["valid_pairs"]), time.perf_counter() - start_time)
    return final, latent_shape, seen, time_values, group_count


def _write_outputs(
    stats: dict[str, Any],
    latent_shape: tuple[int, ...],
    encoded_samples: int,
    time_values: list[Any],
    group_count: dict[tuple[int, ...], int],
    output_dir: Path,
    config: dict[str, Any],
) -> Path:
    local_diag_stats = compute_diag_approx_diagnostics(stats["local_covariance"], stats["local_correlation"])
    global_diag_stats = compute_sampled_diag_approx_diagnostics(
        stats["variance"],
        stats["sampled_offdiag_covariance"],
        stats["sampled_offdiag_correlation"],
        latent_dim=int(np.prod(latent_shape)),
    )
    corr_distribution = compute_offdiag_correlation_distribution(
        stats["sampled_offdiag_correlation"],
        corr_bins=int(config.get("corr_bins", 201)),
        mode="global_all_nmc_pairs" if bool(stats.get("offdiag_pairs_are_exact", False)) else "global_random_nmc_pairs",
    )
    npz_path = output_dir / "latent_nmc_background_covariance.npz"
    LOGGER.info("Saving climatological Bz arrays to %s", npz_path)
    np.savez_compressed(
        npz_path,
        variance=stats["variance"],
        B_diag=stats["B_diag"],
        background_covariance_diag=stats["B_diag"],
        anomaly_mean=stats["perturbation_mean"],
        perturbation_mean=stats["perturbation_mean"],
        offdiag_corr_hist=corr_distribution["offdiag_corr_hist"],
        offdiag_corr_bin_edges=corr_distribution["offdiag_corr_bin_edges"],
        offdiag_abs_corr_hist=corr_distribution["offdiag_abs_corr_hist"],
        offdiag_abs_corr_bin_edges=corr_distribution["offdiag_abs_corr_bin_edges"],
        offdiag_corr_p95=np.array([corr_distribution["offdiag_corr_p95"]], dtype=np.float64),
        offdiag_abs_corr_p95=np.array([corr_distribution["offdiag_abs_corr_p95"]], dtype=np.float64),
        offdiag_corr_pair_count=np.array([corr_distribution["offdiag_corr_pair_count"]], dtype=np.int64),
        offdiag_corr_distribution_mode=np.array([corr_distribution["offdiag_corr_distribution_mode"]], dtype=object),
        sampled_diag_energy=np.array([global_diag_stats["diag_energy"]], dtype=np.float64),
        sampled_offdiag_energy_estimate=np.array([global_diag_stats["offdiag_energy_estimate"]], dtype=np.float64),
        sampled_diag_energy_fraction_estimate=np.array([global_diag_stats["diag_energy_fraction_estimate"]], dtype=np.float64),
        sampled_offdiag_energy_fraction_estimate=np.array([global_diag_stats["offdiag_energy_fraction_estimate"]], dtype=np.float64),
        sampled_offdiag_abs_corr_median=np.array([global_diag_stats["offdiag_abs_corr_median"]], dtype=np.float64),
        sampled_offdiag_abs_corr_p95=np.array([global_diag_stats["offdiag_abs_corr_p95"]], dtype=np.float64),
        local_diag_energy_fraction=np.array([local_diag_stats["diag_energy_fraction"]], dtype=np.float64),
        local_offdiag_energy_fraction=np.array([local_diag_stats["offdiag_energy_fraction"]], dtype=np.float64),
        local_offdiag_abs_corr_median=np.array([local_diag_stats["offdiag_abs_corr_median"]], dtype=np.float64),
        local_offdiag_abs_corr_p95=np.array([local_diag_stats["offdiag_abs_corr_p95"]], dtype=np.float64),
        latent_shape=np.array(latent_shape, dtype=np.int32),
        method=np.array(["climatological_latent_anomaly_covariance"], dtype=object),
        climatology_group=np.array([str(config.get("climatology_group", "month"))], dtype=object),
        climatology_keys=np.array([str(key) for key in sorted(group_count)], dtype=object),
        climatology_group_counts=np.array([group_count[key] for key in sorted(group_count)], dtype=np.int32),
        alpha=np.array([float(config.get("alpha", 1.0))], dtype=np.float64),
        lag=np.array([0], dtype=np.int32),
        split=np.array([str(config.get("split", "all")).lower()], dtype=object),
        val_fraction=np.array([float(config.get("val_fraction", 0.1))], dtype=np.float64),
        diagnostic_pair_count=np.array([corr_distribution["offdiag_corr_pair_count"]], dtype=np.int64),
        diagnostic_seed=np.array([int(config.get("diagnostic_seed", 20260525))], dtype=np.int64),
        diagnostic_pairs_are_exact=np.array([bool(stats.get("offdiag_pairs_are_exact", False))], dtype=bool),
        diagonal_band_radius=np.array([int(config.get("diagonal_band_radius", 0))], dtype=np.int32),
        diagonal_block_size=np.array([int(config.get("diagonal_block_size", config.get("plot_dims", 1)))], dtype=np.int32),
        diagonal_block_starts=stats["diagonal_block_starts"],
        channel_correlation=stats["channel_correlation"],
        spatial_correlation=stats["spatial_correlation"],
        spatial_correlations=stats["spatial_correlations"],
        spatial_channel_index=np.array([stats["spatial_channel_index"]], dtype=np.int32),
        sample_count=np.array([encoded_samples], dtype=np.int32),
        anomaly_count=np.array([stats["valid_pairs"]], dtype=np.int32),
        valid_pairs=np.array([stats["valid_pairs"]], dtype=np.int32),
        time_values=np.array(time_values[:encoded_samples], dtype=object),
    )
    if bool(config.get("save_local_plot_block", False)):
        with np.load(npz_path, allow_pickle=True) as existing:
            payload = {key: existing[key] for key in existing.files}
        payload.update(
            {
                "local_covariance": stats["local_covariance"],
                "local_correlation": stats["local_correlation"],
                "local_eigenvalues": stats["local_eigenvalues"],
                "local_indices": stats["local_indices"],
            }
        )
        np.savez_compressed(npz_path, **payload)
    if bool(config.get("save_diagonal_blocks", False)):
        with np.load(npz_path, allow_pickle=True) as existing:
            payload = {key: existing[key] for key in existing.files}
        payload.update(
            {
                "diagonal_block_covariance": stats["diagonal_block_covariance"],
                "diagonal_block_correlation": stats["diagonal_block_correlation"],
            }
        )
        np.savez_compressed(npz_path, **payload)

    LOGGER.info("Writing climatological Bz figures")
    plot_bdiag(stats["B_diag"], latent_shape, output_dir / "latent_nmc_Bz_diag.png", cmap=str(config.get("diag_cmap", "viridis")))
    if bool(config.get("plot_local_block", False)):
        plot_local_covariance_block(
            stats["local_covariance"],
            stats["local_indices"],
            output_dir / "latent_nmc_Bz_local_block.png",
            title="Local covariance block of climatological $B_z$",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_matrix(
            stats["local_correlation"],
            output_dir / "latent_nmc_local_correlation.png",
            "Local Climatological Correlation Block",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    if bool(config.get("plot_channel_spatial_correlation", False)):
        plot_correlation_matrix(
            stats["channel_correlation"],
            output_dir / "latent_nmc_channel_correlation.png",
            "Channel correlation",
            cmap=str(config.get("matrix_cmap", "bwr")),
            xlabel="Latent channel",
            ylabel="Latent channel",
        )
        plot_correlation_matrix(
            stats["spatial_correlation"],
            output_dir / "latent_nmc_spatial_correlation.png",
            f"Spatial correlation (channel {int(stats['spatial_channel_index'])})",
            cmap=str(config.get("matrix_cmap", "bwr")),
            xlabel="Flattened spatial index",
            ylabel="Flattened spatial index",
        )
        plot_paper_correlation_summary(
            stats["channel_correlation"],
            stats["spatial_correlation"],
            output_dir / "latent_nmc_paper_correlation_summary.png",
            spatial_channel_index=int(stats["spatial_channel_index"]),
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        if bool(config.get("plot_all_spatial_channel_correlations", False)):
            spatial_correlations = np.asarray(stats["spatial_correlations"], dtype=np.float32)
            all_spatial_dir = output_dir / "latent_nmc_spatial_correlations_all"
            all_spatial_dir.mkdir(parents=True, exist_ok=True)
            iterator = range(spatial_correlations.shape[0])
            if tqdm is not None:
                iterator = tqdm(iterator, desc="Plotting spatial channel correlations", unit="channel")
            for channel_index in iterator:
                plot_correlation_matrix(
                    spatial_correlations[channel_index],
                    all_spatial_dir / f"latent_nmc_spatial_correlation_channel_{channel_index:03d}.png",
                    f"Spatial correlation (channel {channel_index})",
                    cmap=str(config.get("matrix_cmap", "bwr")),
                    xlabel="Flattened spatial index",
                    ylabel="Flattened spatial index",
                )
    if bool(config.get("plot_diagonal_band", False)):
        tile_width = int(config.get("diagonal_band_tile_width", 1024))
        plot_diagonal_band_mosaic(
            stats["diagonal_band_covariance"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_Bz_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of climatological $B_z$",
            label="Covariance",
            tile_width=tile_width,
            cmap=str(config.get("matrix_cmap", "bwr")),
            symmetric=True,
        )
        plot_diagonal_band_mosaic(
            stats["diagonal_band_correlation"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_corr_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of climatological correlation",
            label="Correlation",
            tile_width=tile_width,
            cmap=str(config.get("matrix_cmap", "bwr")),
            symmetric=True,
        )
    if bool(config.get("plot_diagonal_blocks", False)) and stats["diagonal_block_covariance"].size:
        tile_cols = int(config.get("diagonal_block_tile_cols", 1))
        plot_diagonal_block_mosaic(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            output_dir / "latent_nmc_Bz_diagonal_block_mosaic.png",
            title="Diagonal blocks of climatological $B_z$",
            label="Covariance",
            tile_cols=tile_cols,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_mosaic(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            output_dir / "latent_nmc_corr_diagonal_block_mosaic.png",
            title="Diagonal blocks of climatological correlation",
            label="Correlation",
            tile_cols=tile_cols,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    if bool(config.get("save_diagonal_block_tiles", False)) and stats["diagonal_block_covariance"].size:
        block_dir = output_dir / "diagonal_block_tiles"
        plot_diagonal_block_tiles(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="Bz_covariance",
            title="Climatological $B_z$ diagonal block",
            label="Covariance",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_tiles(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="correlation",
            title="Climatological correlation diagonal block",
            label="Correlation",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    plot_offdiag_correlation_distribution(
        corr_distribution,
        output_dir / "latent_nmc_offdiag_correlation_distribution.png",
        color=str(config.get("hist_color", "#1f77b4")),
        p95_color=str(config.get("p95_color", "#d88700")),
    )
    LOGGER.info(
        "Offdiag |corr| p95 (%s, pairs=%d): %.4f; sampled diag energy fraction estimate=%.4f; local plot diag fraction=%.4f",
        corr_distribution["offdiag_corr_distribution_mode"],
        int(corr_distribution["offdiag_corr_pair_count"]),
        float(corr_distribution["offdiag_abs_corr_p95"]),
        float(global_diag_stats["diag_energy_fraction_estimate"]),
        float(local_diag_stats["diag_energy_fraction"]),
    )
    LOGGER.info("Climatological Bz output complete: %s", npz_path)
    return npz_path


def run(config: dict[str, Any]) -> Path:
    device = resolve_device(config.get("device"))
    output_dir = _resolve_project_path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = _resolve_project_path(config["checkpoint"])
    dataset = _build_dataset(config)
    max_samples_cfg = config.get("max_samples")
    max_samples = len(dataset) if max_samples_cfg is None else min(int(max_samples_cfg), len(dataset))
    LOGGER.info(
        "Dataset ready: path=%s split=%s samples=%d max_samples=%d device=%s output=%s",
        config["data_path"],
        str(config.get("split", "all")).lower(),
        len(dataset),
        max_samples,
        device,
        output_dir,
    )
    model = load_model(checkpoint, device)
    LOGGER.info("Checkpoint loaded: %s", checkpoint)
    stats, latent_shape, encoded_samples, time_values, group_count = encode_and_compute_climatological_bz(
        model,
        dataset,
        device=device,
        batch_size=int(config.get("batch_size", 1)),
        max_samples=max_samples,
        plot_dims=int(config.get("plot_dims", 1)),
        config=config,
    )
    return _write_outputs(stats, latent_shape, encoded_samples, time_values, group_count, output_dir, config)


def _checkpoint_entries(config: dict, cli_checkpoints: list[str] | None, cli_labels: list[str] | None) -> list[tuple[str, str | None]]:
    if cli_checkpoints:
        labels = cli_labels or []
        return [(checkpoint, labels[index] if index < len(labels) else None) for index, checkpoint in enumerate(cli_checkpoints)]
    entries = config.get("checkpoints")
    if isinstance(entries, list) and entries:
        out = []
        for item in entries:
            if isinstance(item, dict):
                path = item.get("path") or item.get("checkpoint")
                if path:
                    out.append((str(path), item.get("label") or item.get("name")))
            else:
                out.append((str(item), None))
        if out:
            return out
    return [(str(config["checkpoint"]), None)]


def run_versions(config: dict, checkpoints: list[str] | None, labels: list[str] | None) -> list[Path]:
    entries = _checkpoint_entries(config, checkpoints, labels)
    if len(entries) == 1:
        run_config = dict(config)
        run_config["checkpoint"] = entries[0][0]
        return [run(run_config)]
    base_output_dir = _resolve_project_path(config["output_dir"])
    outputs: list[Path] = []
    for checkpoint, label in entries:
        run_config = dict(config)
        version_name = label or _safe_name(checkpoint)
        run_config["checkpoint"] = checkpoint
        run_config["output_dir"] = str(base_output_dir / version_name)
        LOGGER.info("Running climatological Bz version %s -> %s", version_name, run_config["output_dir"])
        outputs.append(run(run_config))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate latent-space Bz from climatology-removed CM2 latent anomalies.")
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "climatological_background_paper.json"))
    parser.add_argument("--checkpoint", nargs="+", default=None, help="Override config checkpoint. Pass multiple paths to write one subdir per version.")
    parser.add_argument("--labels", nargs="+", default=None, help="Optional labels for multiple --checkpoint paths.")
    parser.add_argument("--output-dir", default=None, help="Override config output_dir.")
    parser.add_argument("--data-path", default=None, help="Override config data_path.")
    parser.add_argument("--max-samples", type=int, default=None, help="Override config max_samples.")
    parser.add_argument("--climatology-group", choices=["month", "month_hour", "all"], default=None)
    parser.add_argument("--alpha", type=float, default=None)
    args = parser.parse_args()
    configure_logging()
    with open(args.config, "r", encoding="utf-8") as f:
        config = json.load(f)
    if args.output_dir:
        config["output_dir"] = args.output_dir
    if args.data_path:
        config["data_path"] = args.data_path
    if args.max_samples is not None:
        config["max_samples"] = args.max_samples
    if args.climatology_group:
        config["climatology_group"] = args.climatology_group
    if args.alpha is not None:
        config["alpha"] = args.alpha
    outputs = run_versions(config, args.checkpoint, args.labels)
    print("saved " + ", ".join(str(path) for path in outputs), flush=True)


if __name__ == "__main__":
    main()
