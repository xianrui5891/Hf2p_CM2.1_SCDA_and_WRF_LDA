from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from data_preprocess.climatological_background import (  # noqa: E402
    _build_dataset,
    _climatology_key,
    _resolve_project_path,
)
from data_preprocess.nmc_background import (  # noqa: E402
    OnlineNMCStats,
    _optional_int_list,
    build_loader,
    compute_diag_approx_diagnostics,
    compute_offdiag_correlation_distribution,
    compute_sampled_diag_approx_diagnostics,
    configure_logging,
    plot_diagonal_band_mosaic,
    plot_diagonal_block_mosaic,
    plot_diagonal_block_tiles,
    plot_correlation_matrix,
    plot_local_covariance_block,
    plot_matrix,
    plot_offdiag_correlation_distribution,
    plot_paper_correlation_summary,
    resolve_device,
)

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None


LOGGER = logging.getLogger("cm2_lda.state_climatological_background")
PAPER_DPI = 300


def _time_values(dataset, sample_count: int) -> list[Any]:
    values = list(getattr(dataset, "time_values", list(range(len(dataset)))))
    return values[:sample_count]


def _dataset_vars(dataset) -> tuple[list[str], list[str]]:
    atm_vars = list(getattr(dataset, "atm_vars", [])) or ["TEMP", "UCOMP", "VCOMP", "PS"]
    ocn_vars = list(getattr(dataset, "ocn_vars", [])) or ["SST", "U_SURF", "V_SURF", "ETA_T"]
    return atm_vars, ocn_vars


def _unpack_state_batch(batch) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None, torch.Tensor | None]:
    if not isinstance(batch, (tuple, list)) or len(batch) < 2:
        raise ValueError(f"State dataset batch must contain at least atm and ocn tensors, got {type(batch).__name__}.")
    atm = batch[0]
    ocn = batch[1]
    atm_mask = batch[2] if len(batch) > 2 else None
    ocn_mask = batch[3] if len(batch) > 3 else None
    return atm, ocn, atm_mask, ocn_mask


def _masked_to_device(
    values: torch.Tensor,
    mask: torch.Tensor | None,
    keep: int,
    device: torch.device,
) -> torch.Tensor:
    out = torch.nan_to_num(values[:keep].to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
    if mask is not None:
        valid = mask[:keep].to(device, non_blocking=True).bool()
        out = torch.where(valid, out, torch.zeros_like(out))
    return out


def _flatten_state_batch(batch, keep: int, device: torch.device) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...]]:
    atm, ocn, atm_mask, ocn_mask = _unpack_state_batch(batch)
    atm = _masked_to_device(atm, atm_mask, keep, device)
    ocn = _masked_to_device(ocn, ocn_mask, keep, device)
    atm_shape = tuple(int(item) for item in atm.shape[1:])
    ocn_shape = tuple(int(item) for item in ocn.shape[1:])
    flat = torch.cat((atm.reshape(atm.shape[0], -1), ocn.reshape(ocn.shape[0], -1)), dim=1)
    return flat.float().cpu().numpy(), atm_shape, ocn_shape


def _make_stats(state_dim: int, plot_dims: int, config: dict[str, Any]) -> OnlineNMCStats:
    return OnlineNMCStats(
        latent_dim=state_dim,
        latent_shape=None,
        local_dims=plot_dims,
        diagnostic_pair_count=int(config.get("diagnostic_pair_count", 1_000_000)),
        diagnostic_seed=int(config.get("diagnostic_seed", 20260525)),
        pair_chunk_size=int(config.get("diagnostic_pair_chunk_size", 200_000)),
        diagonal_band_radius=int(config.get("diagonal_band_radius", 0)),
        diagonal_block_size=int(config.get("diagonal_block_size", config.get("plot_dims", 1))),
        diagonal_block_count=int(config.get("diagonal_block_count", 0)),
        diagonal_block_starts=_optional_int_list(config.get("diagonal_block_starts")),
        channel_spatial_diagnostics=False,
        spatial_channel_index=0,
        all_spatial_channel_correlations=False,
    )


def _correlation_from_sums(sum_: np.ndarray, outer: np.ndarray, count: int) -> np.ndarray:
    if count < 2 or outer.size == 0:
        return np.empty_like(outer, dtype=np.float32)
    mean = sum_ / count
    covariance = (outer - count * np.outer(mean, mean)) / (count - 1)
    covariance = np.nan_to_num(covariance, nan=0.0, posinf=0.0, neginf=0.0)
    std = np.sqrt(np.clip(np.diag(covariance), 1e-24, None))
    correlation = covariance / np.outer(std, std)
    correlation = np.clip(np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
    return correlation.astype(np.float32)


def _parse_spatial_variable(value: str, atm_vars: list[str], ocn_vars: list[str]) -> tuple[str, int, str]:
    text = str(value or "ocn:SST")
    if ":" in text:
        domain, name = text.split(":", 1)
        domain = domain.strip().lower()
        name = name.strip()
    else:
        name = text.strip()
        domain = "ocn" if name in ocn_vars else "atm"
    vars_ = ocn_vars if domain == "ocn" else atm_vars
    lookup = {item.upper(): idx for idx, item in enumerate(vars_)}
    index = lookup.get(name.upper(), 0)
    return domain, index, vars_[index]


class StateCorrelationDiagnostics:
    def __init__(
        self,
        atm_shape: tuple[int, ...],
        ocn_shape: tuple[int, ...],
        atm_vars: list[str],
        ocn_vars: list[str],
        config: dict[str, Any],
    ) -> None:
        self.enabled = bool(config.get("channel_spatial_diagnostics", False))
        self.all_spatial = bool(config.get("plot_all_spatial_channel_correlations", False))
        self.atm_shape = tuple(int(item) for item in atm_shape)
        self.ocn_shape = tuple(int(item) for item in ocn_shape)
        self.atm_vars = list(atm_vars)
        self.ocn_vars = list(ocn_vars)
        self.variable_names = [f"atm:{name}" for name in atm_vars] + [f"ocn:{name}" for name in ocn_vars]
        self.atm_channels = self.atm_shape[0]
        self.ocn_channels = self.ocn_shape[0]
        self.channel_count = self.atm_channels + self.ocn_channels
        self.channel_correlation = np.empty((0, 0), dtype=np.float32)
        self.spatial_correlation = np.empty((0, 0), dtype=np.float32)
        self.spatial_correlations_by_name: dict[str, np.ndarray] = {}
        self.selected_domain, self.selected_index, self.selected_name = _parse_spatial_variable(
            str(config.get("spatial_variable", "ocn:SST")),
            self.atm_vars,
            self.ocn_vars,
        )
        self.selected_label = f"{self.selected_domain}:{self.selected_name}"
        self.atm_stride = max(int(config.get("atm_spatial_stride", 4)), 1)
        self.ocn_stride = max(int(config.get("ocn_spatial_stride", 8)), 1)
        if not self.enabled:
            return

        self.atm_channel_count = 0
        self.atm_channel_sum = np.zeros(self.atm_channels, dtype=np.float64)
        self.atm_channel_outer = np.zeros((self.atm_channels, self.atm_channels), dtype=np.float64)
        self.ocn_channel_count = 0
        self.ocn_channel_sum = np.zeros(self.ocn_channels, dtype=np.float64)
        self.ocn_channel_outer = np.zeros((self.ocn_channels, self.ocn_channels), dtype=np.float64)

        self.spatial_sum: dict[str, np.ndarray] = {}
        self.spatial_outer: dict[str, np.ndarray] = {}
        self.spatial_count: dict[str, int] = {}

    def _downsample(self, field: np.ndarray, domain: str) -> np.ndarray:
        stride = self.ocn_stride if domain == "ocn" else self.atm_stride
        return field[:, ::stride, ::stride].reshape(field.shape[0], -1)

    def _update_spatial(self, label: str, rows: np.ndarray) -> None:
        rows = np.asarray(rows, dtype=np.float64)
        if label not in self.spatial_sum:
            self.spatial_sum[label] = np.zeros(rows.shape[1], dtype=np.float64)
            self.spatial_outer[label] = np.zeros((rows.shape[1], rows.shape[1]), dtype=np.float64)
            self.spatial_count[label] = 0
        self.spatial_sum[label] += rows.sum(axis=0)
        self.spatial_outer[label] += rows.T @ rows
        self.spatial_count[label] += int(rows.shape[0])

    def update(self, atm_anomaly: np.ndarray, ocn_anomaly: np.ndarray) -> None:
        if not self.enabled:
            return
        atm = np.asarray(atm_anomaly, dtype=np.float64).reshape((-1,) + self.atm_shape)
        ocn = np.asarray(ocn_anomaly, dtype=np.float64).reshape((-1,) + self.ocn_shape)

        atm_rows = np.moveaxis(atm, 1, -1).reshape(-1, self.atm_channels)
        self.atm_channel_sum += atm_rows.sum(axis=0)
        self.atm_channel_outer += atm_rows.T @ atm_rows
        self.atm_channel_count += int(atm_rows.shape[0])

        ocn_rows = np.moveaxis(ocn, 1, -1).reshape(-1, self.ocn_channels)
        self.ocn_channel_sum += ocn_rows.sum(axis=0)
        self.ocn_channel_outer += ocn_rows.T @ ocn_rows
        self.ocn_channel_count += int(ocn_rows.shape[0])

        if self.selected_domain == "atm":
            self._update_spatial(self.selected_label, self._downsample(atm[:, self.selected_index], "atm"))
        else:
            self._update_spatial(self.selected_label, self._downsample(ocn[:, self.selected_index], "ocn"))

        if self.all_spatial:
            for index, name in enumerate(self.atm_vars):
                label = f"atm:{name}"
                if label != self.selected_label:
                    self._update_spatial(label, self._downsample(atm[:, index], "atm"))
            for index, name in enumerate(self.ocn_vars):
                label = f"ocn:{name}"
                if label != self.selected_label:
                    self._update_spatial(label, self._downsample(ocn[:, index], "ocn"))

    def finalize(self) -> dict[str, Any]:
        if not self.enabled:
            return {
                "channel_correlation": np.empty((0, 0), dtype=np.float32),
                "spatial_correlation": np.empty((0, 0), dtype=np.float32),
                "spatial_correlations_by_name": {},
                "selected_spatial_label": self.selected_label,
                "variable_names": self.variable_names,
            }
        atm_corr = _correlation_from_sums(self.atm_channel_sum, self.atm_channel_outer, self.atm_channel_count)
        ocn_corr = _correlation_from_sums(self.ocn_channel_sum, self.ocn_channel_outer, self.ocn_channel_count)
        channel_corr = np.zeros((self.channel_count, self.channel_count), dtype=np.float32)
        channel_corr[: self.atm_channels, : self.atm_channels] = atm_corr
        channel_corr[self.atm_channels :, self.atm_channels :] = ocn_corr

        spatial_by_name = {
            label: _correlation_from_sums(self.spatial_sum[label], self.spatial_outer[label], self.spatial_count[label])
            for label in sorted(self.spatial_sum)
        }
        selected_spatial = spatial_by_name.get(self.selected_label, np.empty((0, 0), dtype=np.float32))
        return {
            "channel_correlation": channel_corr,
            "spatial_correlation": selected_spatial,
            "spatial_correlations_by_name": spatial_by_name,
            "selected_spatial_label": self.selected_label,
            "variable_names": self.variable_names,
        }


def _split_bdiag(
    bdiag: np.ndarray,
    atm_shape: tuple[int, ...],
    ocn_shape: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    atm_dim = int(np.prod(atm_shape))
    atm = np.asarray(bdiag[:atm_dim], dtype=np.float32).reshape(atm_shape)
    ocn = np.asarray(bdiag[atm_dim : atm_dim + int(np.prod(ocn_shape))], dtype=np.float32).reshape(ocn_shape)
    return atm, ocn


def _plot_domain_bdiag(fields: np.ndarray, var_names: list[str], output_path: Path, title: str, cmap: str) -> None:
    channel_count = int(fields.shape[0])
    cols = min(channel_count, 4)
    rows = int(math.ceil(channel_count / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.4 * rows), squeeze=False, constrained_layout=True)
    finite = fields[np.isfinite(fields)]
    vmax = float(np.percentile(finite, 99.5)) if finite.size else 1.0
    vmax = max(vmax, 1e-12)
    image = None
    for idx in range(rows * cols):
        ax = axes[idx // cols, idx % cols]
        if idx >= channel_count:
            ax.axis("off")
            continue
        image = ax.imshow(fields[idx], origin="lower", cmap=cmap, vmin=0.0, vmax=vmax)
        label = var_names[idx] if idx < len(var_names) else f"var {idx}"
        ax.set_title(label, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
    if image is not None:
        fig.colorbar(image, ax=axes, shrink=0.82, label="B diagonal variance")
    fig.suptitle(title, fontsize=13)
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def _plot_combined_bdiag(
    bdiag: np.ndarray,
    atm_shape: tuple[int, ...],
    ocn_shape: tuple[int, ...],
    atm_vars: list[str],
    ocn_vars: list[str],
    output_path: Path,
    cmap: str,
) -> None:
    atm_fields, ocn_fields = _split_bdiag(bdiag, atm_shape, ocn_shape)
    atm_mean = atm_fields.mean(axis=0)
    ocn_mean = ocn_fields.mean(axis=0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, field, title in (
        (axes[0], atm_mean, "Atmosphere mean B diag"),
        (axes[1], ocn_mean, "Ocean mean B diag"),
    ):
        finite = field[np.isfinite(field)]
        vmax = float(np.percentile(finite, 99.5)) if finite.size else 1.0
        image = ax.imshow(field, origin="lower", cmap=cmap, vmin=0.0, vmax=max(vmax, 1e-12))
        ax.set_title(title, fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, shrink=0.82)
    fig.suptitle("Direct variable-space climatological $B$ diagonal", fontsize=14)
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)
    _plot_domain_bdiag(atm_fields, atm_vars, output_path.with_name("state_nmc_Bz_atm_variables.png"), "Atmosphere variable-space B diagonal", cmap)
    _plot_domain_bdiag(ocn_fields, ocn_vars, output_path.with_name("state_nmc_Bz_ocn_variables.png"), "Ocean variable-space B diagonal", cmap)


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
            arr = np.asarray(stats[key])
            stats[key] = (arr * alpha2).astype(arr.dtype, copy=False)
    return stats


@torch.inference_mode()
def encode_and_compute_state_b(
    dataset,
    device: torch.device,
    batch_size: int,
    max_samples: int,
    plot_dims: int,
    config: dict[str, Any],
) -> tuple[dict[str, Any], tuple[int, ...], tuple[int, ...], int, list[Any], dict[tuple[int, ...], int]]:
    total_samples = min(len(dataset), max_samples)
    time_values = _time_values(dataset, total_samples)
    loader = build_loader(dataset, batch_size=batch_size, config=config)
    group_sum: dict[tuple[int, ...], np.ndarray] = {}
    group_count: dict[tuple[int, ...], int] = {}
    atm_shape: tuple[int, ...] | None = None
    ocn_shape: tuple[int, ...] | None = None
    state_dim = 0
    seen = 0
    start_time = time.perf_counter()
    LOGGER.info(
        "State climatological B start: samples=%d group=%s batch_size=%d loader_workers=%d",
        total_samples,
        config.get("climatology_group", "month_hour"),
        batch_size,
        int(config.get("num_workers", min(4, os.cpu_count() or 4))),
    )
    progress = tqdm(total=total_samples, desc="Pass 1/2 state climatology", unit="sample") if tqdm is not None else None
    for batch in loader:
        if seen >= total_samples:
            break
        atm, _, _, _ = _unpack_state_batch(batch)
        keep = min(int(atm.shape[0]), total_samples - seen)
        if keep <= 0:
            break
        flat, current_atm_shape, current_ocn_shape = _flatten_state_batch(batch, keep, device)
        if atm_shape is None:
            atm_shape = current_atm_shape
            ocn_shape = current_ocn_shape
            state_dim = int(flat.shape[1])
            LOGGER.info("State shape: atm=%s ocn=%s state_dim=%d", atm_shape, ocn_shape, state_dim)
        for local_index, row in enumerate(flat):
            sample_index = seen + local_index
            key = _climatology_key(time_values[sample_index], sample_index, config)
            if key not in group_sum:
                group_sum[key] = np.zeros(state_dim, dtype=np.float64)
                group_count[key] = 0
            group_sum[key] += row.astype(np.float64, copy=False)
            group_count[key] += 1
        seen += keep
        if progress is not None:
            progress.update(keep)
    if progress is not None:
        progress.close()
    if atm_shape is None or ocn_shape is None:
        raise RuntimeError("No samples were read for state climatology statistics.")
    climatology = {
        key: (group_sum[key] / max(group_count[key], 1)).astype(np.float32)
        for key in sorted(group_sum)
    }
    LOGGER.info("State climatology groups: %s", {str(key): group_count[key] for key in sorted(group_count)})

    loader = build_loader(dataset, batch_size=batch_size, config=config)
    stats = _make_stats(state_dim, plot_dims, config)
    diagnostics = StateCorrelationDiagnostics(atm_shape, ocn_shape, *_dataset_vars(dataset), config)
    seen = 0
    log_every_batches = int(config.get("log_every_n_batches", 200))
    progress = tqdm(total=total_samples, desc="Pass 2/2 state anomalies", unit="sample") if tqdm is not None else None
    for batch_index, batch in enumerate(loader, start=1):
        if seen >= total_samples:
            break
        atm, _, _, _ = _unpack_state_batch(batch)
        keep = min(int(atm.shape[0]), total_samples - seen)
        if keep <= 0:
            break
        flat, _, _ = _flatten_state_batch(batch, keep, device)
        anomalies = np.empty_like(flat, dtype=np.float32)
        for local_index, row in enumerate(flat):
            sample_index = seen + local_index
            key = _climatology_key(time_values[sample_index], sample_index, config)
            anomalies[local_index] = row - climatology[key]
        stats.update(anomalies)
        atm_dim = int(np.prod(atm_shape))
        diagnostics.update(
            anomalies[:, :atm_dim].reshape((anomalies.shape[0],) + atm_shape),
            anomalies[:, atm_dim:].reshape((anomalies.shape[0],) + ocn_shape),
        )
        seen += keep
        if progress is not None:
            progress.update(keep)
        if batch_index % log_every_batches == 0:
            elapsed = max(time.perf_counter() - start_time, 1e-6)
            LOGGER.info("State anomalies processed %d/%d, rate=%.2f sample/s", seen, total_samples, seen / elapsed)
    if progress is not None:
        progress.close()
    final = stats.finalize()
    final.update(diagnostics.finalize())
    final["sample_count"] = int(seen)
    final = _scale_covariance_outputs(final, float(config.get("alpha", 1.0)))
    LOGGER.info("State climatological B done: samples=%d anomalies=%d elapsed=%.1fs", seen, int(final["valid_pairs"]), time.perf_counter() - start_time)
    return final, atm_shape, ocn_shape, seen, time_values, group_count


def _write_outputs(
    stats: dict[str, Any],
    atm_shape: tuple[int, ...],
    ocn_shape: tuple[int, ...],
    atm_vars: list[str],
    ocn_vars: list[str],
    encoded_samples: int,
    time_values: list[Any],
    group_count: dict[tuple[int, ...], int],
    output_dir: Path,
    config: dict[str, Any],
) -> Path:
    state_shape = np.array([int(np.prod(atm_shape)) + int(np.prod(ocn_shape))], dtype=np.int32)
    local_diag_stats = compute_diag_approx_diagnostics(stats["local_covariance"], stats["local_correlation"])
    global_diag_stats = compute_sampled_diag_approx_diagnostics(
        stats["variance"],
        stats["sampled_offdiag_covariance"],
        stats["sampled_offdiag_correlation"],
        latent_dim=int(state_shape[0]),
    )
    corr_distribution = compute_offdiag_correlation_distribution(
        stats["sampled_offdiag_correlation"],
        corr_bins=int(config.get("corr_bins", 201)),
        mode="global_all_nmc_pairs" if bool(stats.get("offdiag_pairs_are_exact", False)) else "global_random_nmc_pairs",
    )
    npz_path = output_dir / "latent_nmc_background_covariance.npz"
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
        latent_shape=state_shape,
        state_shape=state_shape,
        atm_shape=np.array(atm_shape, dtype=np.int32),
        ocn_shape=np.array(ocn_shape, dtype=np.int32),
        atm_dim=np.array([int(np.prod(atm_shape))], dtype=np.int32),
        ocn_dim=np.array([int(np.prod(ocn_shape))], dtype=np.int32),
        atm_vars=np.array(atm_vars, dtype=object),
        ocn_vars=np.array(ocn_vars, dtype=object),
        method=np.array(["direct_variable_climatological_anomaly_covariance"], dtype=object),
        climatology_group=np.array([str(config.get("climatology_group", "month_hour"))], dtype=object),
        climatology_keys=np.array([str(key) for key in sorted(group_count)], dtype=object),
        climatology_group_counts=np.array([group_count[key] for key in sorted(group_count)], dtype=np.int32),
        alpha=np.array([float(config.get("alpha", 1.0))], dtype=np.float64),
        lag=np.array([0], dtype=np.int32),
        split=np.array([str(config.get("split", "val")).lower()], dtype=object),
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
        selected_spatial_label=np.array([stats.get("selected_spatial_label", "")], dtype=object),
        state_variable_names=np.array(stats.get("variable_names", []), dtype=object),
        spatial_channel_index=np.array([0], dtype=np.int32),
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

    cmap = str(config.get("diag_cmap", "viridis"))
    _plot_combined_bdiag(
        stats["B_diag"],
        atm_shape,
        ocn_shape,
        atm_vars,
        ocn_vars,
        output_dir / "latent_nmc_Bz_diag.png",
        cmap,
    )
    if bool(config.get("plot_channel_spatial_correlation", False)):
        plot_correlation_matrix(
            stats["channel_correlation"],
            output_dir / "latent_nmc_channel_correlation.png",
            "Variable channel correlation",
            cmap=str(config.get("matrix_cmap", "bwr")),
            xlabel="Variable channel",
            ylabel="Variable channel",
        )
        plot_correlation_matrix(
            stats["spatial_correlation"],
            output_dir / "latent_nmc_spatial_correlation.png",
            f"Spatial correlation ({stats.get('selected_spatial_label', 'selected')})",
            cmap=str(config.get("matrix_cmap", "bwr")),
            xlabel="Downsampled spatial index",
            ylabel="Downsampled spatial index",
        )
        plot_paper_correlation_summary(
            stats["channel_correlation"],
            stats["spatial_correlation"],
            output_dir / "latent_nmc_paper_correlation_summary.png",
            spatial_channel_index=0,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        if bool(config.get("plot_all_spatial_channel_correlations", False)):
            all_spatial_dir = output_dir / "latent_nmc_spatial_correlations_all"
            all_spatial_dir.mkdir(parents=True, exist_ok=True)
            spatial_by_name = dict(stats.get("spatial_correlations_by_name", {}))
            iterator = sorted(spatial_by_name)
            if tqdm is not None:
                iterator = tqdm(iterator, desc="Plotting variable spatial correlations", unit="variable")
            for label in iterator:
                safe_label = label.replace(":", "_").replace("/", "_")
                plot_correlation_matrix(
                    spatial_by_name[label],
                    all_spatial_dir / f"latent_nmc_spatial_correlation_{safe_label}.png",
                    f"Spatial correlation ({label})",
                    cmap=str(config.get("matrix_cmap", "bwr")),
                    xlabel="Downsampled spatial index",
                    ylabel="Downsampled spatial index",
                )
    if bool(config.get("plot_local_block", False)):
        plot_local_covariance_block(
            stats["local_covariance"],
            stats["local_indices"],
            output_dir / "latent_nmc_Bz_local_block.png",
            title="Local covariance block of direct variable-space $B$",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_matrix(
            stats["local_correlation"],
            output_dir / "latent_nmc_local_correlation.png",
            "Local Direct Variable-Space Correlation Block",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    if bool(config.get("plot_diagonal_band", False)):
        tile_width = int(config.get("diagonal_band_tile_width", 1024))
        plot_diagonal_band_mosaic(
            stats["diagonal_band_covariance"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_Bz_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of direct variable-space $B$",
            label="Covariance",
            tile_width=tile_width,
            cmap=str(config.get("matrix_cmap", "bwr")),
            symmetric=True,
        )
        plot_diagonal_band_mosaic(
            stats["diagonal_band_correlation"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_corr_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of direct variable-space correlation",
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
            title="Diagonal blocks of direct variable-space $B$",
            label="Covariance",
            tile_cols=tile_cols,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_mosaic(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            output_dir / "latent_nmc_corr_diagonal_block_mosaic.png",
            title="Diagonal blocks of direct variable-space correlation",
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
            title="Direct variable-space $B$ diagonal block",
            label="Covariance",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_tiles(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="correlation",
            title="Direct variable-space correlation diagonal block",
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
        "Direct variable-space output complete: %s; offdiag |corr| p95 (%s, pairs=%d)=%.4f",
        npz_path,
        corr_distribution["offdiag_corr_distribution_mode"],
        int(corr_distribution["offdiag_corr_pair_count"]),
        float(corr_distribution["offdiag_abs_corr_p95"]),
    )
    return npz_path


def run(config: dict[str, Any]) -> Path:
    device = resolve_device(config.get("device"))
    output_dir = _resolve_project_path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset = _build_dataset(config)
    atm_vars, ocn_vars = _dataset_vars(dataset)
    max_samples_cfg = config.get("max_samples")
    max_samples = len(dataset) if max_samples_cfg is None else min(int(max_samples_cfg), len(dataset))
    LOGGER.info(
        "Dataset ready: path=%s split=%s samples=%d max_samples=%d vars=%s/%s device=%s output=%s",
        config["data_path"],
        str(config.get("split", "val")).lower(),
        len(dataset),
        max_samples,
        atm_vars,
        ocn_vars,
        device,
        output_dir,
    )
    stats, atm_shape, ocn_shape, encoded_samples, time_values, group_count = encode_and_compute_state_b(
        dataset,
        device=device,
        batch_size=int(config.get("batch_size", 1)),
        max_samples=max_samples,
        plot_dims=int(config.get("plot_dims", 1)),
        config=config,
    )
    return _write_outputs(
        stats,
        atm_shape,
        ocn_shape,
        atm_vars,
        ocn_vars,
        encoded_samples,
        time_values,
        group_count,
        output_dir,
        config,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate direct CM2 variable-space B from climatology-removed anomalies.")
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "state_climatological_background_paper.json"))
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
    output = run(config)
    print(f"saved {output}", flush=True)


if __name__ == "__main__":
    main()
