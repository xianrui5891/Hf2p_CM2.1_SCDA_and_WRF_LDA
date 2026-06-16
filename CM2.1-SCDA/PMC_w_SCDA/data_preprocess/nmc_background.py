from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from collections import deque
from functools import partial
from inspect import signature
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import torch
from torch.utils.data import DataLoader

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from model import CoupledAE
from utils.data import CM2DatasetView, CM2FrameDataset, CM2TensorDataset, FRAME_MANIFEST

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - tqdm is optional for bare environments.
    tqdm = None


LOGGER = logging.getLogger("cm2_lda.nmc_background")
DEFAULT_MAX_SAMPLES = 10_000
_DATALOADER_PARAMS = set(signature(DataLoader).parameters)
PAPER_DPI = 300


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )


def resolve_device(requested: str | None) -> torch.device:
    if requested:
        return torch.device(requested)
    if not torch.cuda.is_available():
        return torch.device("cpu")
    best_index = 0
    best_free_mem = -1
    for index in range(torch.cuda.device_count()):
        try:
            with torch.cuda.device(index):
                free_mem, _ = torch.cuda.mem_get_info()
        except Exception:
            free_mem = -1
        if free_mem > best_free_mem:
            best_free_mem = free_mem
            best_index = index
    return torch.device(f"cuda:{best_index}")


def load_model(checkpoint_path: str | Path, device: torch.device) -> CoupledAE:
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model_cfg = ckpt.get("model_config", {})
    width = int(model_cfg.get("width", 64))
    model = CoupledAE(
        latent_dim=model_cfg.get("latent_dim", ckpt.get("latent_dim")),
        latent_channels=model_cfg.get("latent_channels"),
        in_atm_channels=model_cfg.get("in_atm_channels", 4),
        in_ocn_channels=model_cfg.get("in_ocn_channels", 4),
        width=width,
        decoder_width=model_cfg.get("decoder_width", width * 2),
        decoder_middle_depth=model_cfg.get("decoder_middle_depth", 2),
        decoder_stage_depth=model_cfg.get("decoder_stage_depth", 1),
        atm_decoder_stage_depth=model_cfg.get("atm_decoder_stage_depth"),
        ocn_decoder_stage_depth=model_cfg.get("ocn_decoder_stage_depth"),
        cross_depth=model_cfg.get("cross_depth", 4),
        heads=model_cfg.get("heads", 4),
        attention_backend=model_cfg.get("attention_backend", "auto"),
        fusion_mode=model_cfg.get("fusion_mode", "cross_attention"),
        cnn_fusion_depth=int(model_cfg.get("cnn_fusion_depth", 2)),
    ).to(device)
    state = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt))
    state = {key.removeprefix("module."): value for key, value in state.items()}
    incompatible = model.load_state_dict(state, strict=False)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        LOGGER.warning(
            "Loaded checkpoint with non-strict state dict: missing=%d unexpected=%d",
            len(incompatible.missing_keys),
            len(incompatible.unexpected_keys),
        )
    model.eval()
    return model


def _unpack_atm_ocn_batch(batch) -> tuple[torch.Tensor, torch.Tensor]:
    if not isinstance(batch, (tuple, list)) or len(batch) < 2:
        raise ValueError(f"NMC dataset batch must contain at least atm and ocn tensors, got {type(batch).__name__}.")
    return batch[0], batch[1]


def _set_worker_torch_threads(num_threads: int, worker_id: int) -> None:
    del worker_id
    if num_threads > 0:
        torch.set_num_threads(num_threads)


def build_loader(dataset, batch_size: int, config: dict[str, Any]) -> DataLoader:
    num_workers = int(config.get("num_workers", min(4, os.cpu_count() or 4)))
    pin_memory = bool(config.get("pin_memory", torch.cuda.is_available()))
    kwargs: dict[str, Any] = {
        "batch_size": batch_size,
        "shuffle": False,
        "num_workers": num_workers,
        "pin_memory": pin_memory,
    }
    if num_workers > 0:
        kwargs["persistent_workers"] = bool(config.get("persistent_workers", True))
        kwargs["prefetch_factor"] = int(config.get("prefetch_factor", 8))
        worker_torch_threads = int(config.get("worker_torch_threads", 1))
        kwargs["worker_init_fn"] = partial(_set_worker_torch_threads, worker_torch_threads)
        # NMC lag pairs require deterministic sample order, so keep in-order loading.
        if "in_order" in _DATALOADER_PARAMS:
            kwargs["in_order"] = bool(config.get("in_order", True))
    timeout = int(config.get("timeout", 0))
    if timeout > 0:
        kwargs["timeout"] = timeout
    return DataLoader(dataset, **kwargs)


def _optional_int_list(value: Any) -> list[int] | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [int(item) for item in value]
    return [int(value)]


def _autocast_context(device: torch.device, config: dict[str, Any]):
    amp_enabled = bool(config.get("amp_enabled", device.type == "cuda"))
    dtype_name = str(config.get("amp_dtype", "float16")).lower()
    dtype = torch.bfloat16 if dtype_name in {"bf16", "bfloat16"} else torch.float16
    return torch.amp.autocast(device_type=device.type, dtype=dtype, enabled=amp_enabled and device.type == "cuda")


def _make_offdiag_pairs(latent_dim: int, pair_count: int, seed: int) -> tuple[np.ndarray, np.ndarray, bool]:
    if latent_dim < 2 or pair_count <= 0:
        return np.empty((0,), dtype=np.int64), np.empty((0,), dtype=np.int64), False
    full_pair_count = int(latent_dim) * int(latent_dim - 1)
    if full_pair_count <= int(pair_count):
        row = np.repeat(np.arange(latent_dim, dtype=np.int64), latent_dim - 1)
        col = np.empty(full_pair_count, dtype=np.int64)
        cursor = 0
        base = np.arange(latent_dim, dtype=np.int64)
        for i in range(latent_dim):
            values = base[base != i]
            col[cursor : cursor + values.size] = values
            cursor += values.size
        return row, col, True
    rng = np.random.default_rng(int(seed))
    row = rng.integers(0, latent_dim, size=int(pair_count), dtype=np.int64)
    offset = rng.integers(1, latent_dim, size=int(pair_count), dtype=np.int64)
    col = (row + offset) % latent_dim
    return row, col, False


def _make_diagonal_block_starts(
    latent_dim: int,
    block_size: int,
    block_count: int,
    requested_starts: list[int] | None = None,
) -> np.ndarray:
    latent_dim = int(latent_dim)
    block_size = max(int(block_size), 1)
    if latent_dim <= 0:
        return np.empty((0,), dtype=np.int32)
    max_start = max(latent_dim - block_size, 0)
    if requested_starts:
        starts = sorted({min(max(int(start), 0), max_start) for start in requested_starts})
        return np.asarray(starts, dtype=np.int32)
    block_count = max(int(block_count), 0)
    if block_count == 0:
        return np.empty((0,), dtype=np.int32)
    if block_count == 1 or max_start == 0:
        return np.asarray([0], dtype=np.int32)
    starts = np.linspace(0, max_start, block_count)
    return np.unique(np.rint(starts).astype(np.int32))


def _covariance_from_sums(sum_: np.ndarray, outer: np.ndarray, count: int) -> np.ndarray:
    if count < 2:
        return np.zeros_like(outer, dtype=np.float64)
    mean = sum_ / count
    covariance = (outer - count * np.outer(mean, mean)) / (count - 1)
    return np.nan_to_num(covariance, nan=0.0, posinf=0.0, neginf=0.0)


def _correlation_from_covariance(covariance: np.ndarray) -> np.ndarray:
    if covariance.size == 0:
        return np.empty_like(covariance, dtype=np.float64)
    diag = np.sqrt(np.clip(np.diag(covariance), 1e-12, None))
    correlation = covariance / np.outer(diag, diag)
    return np.clip(np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)


class OnlineNMCStats:
    def __init__(
        self,
        latent_dim: int,
        latent_shape: tuple[int, ...] | None,
        local_dims: int,
        *,
        diagnostic_pair_count: int,
        diagnostic_seed: int,
        pair_chunk_size: int,
        diagonal_band_radius: int,
        diagonal_block_size: int,
        diagonal_block_count: int,
        diagonal_block_starts: list[int] | None,
        channel_spatial_diagnostics: bool,
        spatial_channel_index: int,
        all_spatial_channel_correlations: bool,
    ) -> None:
        self.latent_dim = int(latent_dim)
        self.latent_shape = tuple(int(item) for item in latent_shape) if latent_shape is not None else None
        self.channel_spatial_diagnostics = bool(channel_spatial_diagnostics and self.latent_shape and len(self.latent_shape) == 3)
        self.all_spatial_channel_correlations = bool(all_spatial_channel_correlations and self.channel_spatial_diagnostics)
        if self.channel_spatial_diagnostics:
            channels, height, width = self.latent_shape
            self.channel_count = int(channels)
            self.spatial_dim = int(height * width)
            self.spatial_channel_index = int(np.clip(int(spatial_channel_index), 0, self.channel_count - 1))
            self.channel_sample_count = 0
            self.channel_sum = np.zeros(self.channel_count, dtype=np.float64)
            self.channel_outer = np.zeros((self.channel_count, self.channel_count), dtype=np.float64)
            self.spatial_sample_count = 0
            self.spatial_sum = np.zeros(self.spatial_dim, dtype=np.float64)
            self.spatial_outer = np.zeros((self.spatial_dim, self.spatial_dim), dtype=np.float64)
            self.all_spatial_sample_count = 0
            if self.all_spatial_channel_correlations:
                self.all_spatial_sum = np.zeros((self.channel_count, self.spatial_dim), dtype=np.float64)
                self.all_spatial_outer = np.zeros(
                    (self.channel_count, self.spatial_dim, self.spatial_dim),
                    dtype=np.float64,
                )
            else:
                self.all_spatial_sum = np.empty((0, 0), dtype=np.float64)
                self.all_spatial_outer = np.empty((0, 0, 0), dtype=np.float64)
        else:
            self.channel_count = 0
            self.spatial_dim = 0
            self.spatial_channel_index = 0
            self.channel_sample_count = 0
            self.channel_sum = np.empty((0,), dtype=np.float64)
            self.channel_outer = np.empty((0, 0), dtype=np.float64)
            self.spatial_sample_count = 0
            self.spatial_sum = np.empty((0,), dtype=np.float64)
            self.spatial_outer = np.empty((0, 0), dtype=np.float64)
            self.all_spatial_sample_count = 0
            self.all_spatial_sum = np.empty((0, 0), dtype=np.float64)
            self.all_spatial_outer = np.empty((0, 0, 0), dtype=np.float64)
        self.local_dims = int(min(local_dims, latent_dim))
        self.pair_chunk_size = max(int(pair_chunk_size), 1)
        self.diagonal_band_radius = max(int(diagonal_band_radius), 0)
        self.count = 0
        self.sum = np.zeros(self.latent_dim, dtype=np.float64)
        self.sumsq = np.zeros(self.latent_dim, dtype=np.float64)
        self.local_sum = np.zeros(self.local_dims, dtype=np.float64)
        self.local_outer = np.zeros((self.local_dims, self.local_dims), dtype=np.float64)
        self.diagonal_band_cross_sum = np.zeros(
            (self.diagonal_band_radius + 1, self.latent_dim),
            dtype=np.float64,
        )
        self.diagonal_block_size = min(max(int(diagonal_block_size), 1), self.latent_dim)
        self.diagonal_block_starts = _make_diagonal_block_starts(
            self.latent_dim,
            self.diagonal_block_size,
            int(diagonal_block_count),
            diagonal_block_starts,
        )
        self.diagonal_block_sum = np.zeros(
            (self.diagonal_block_starts.size, self.diagonal_block_size),
            dtype=np.float64,
        )
        self.diagonal_block_outer = np.zeros(
            (self.diagonal_block_starts.size, self.diagonal_block_size, self.diagonal_block_size),
            dtype=np.float64,
        )
        self.offdiag_pair_i, self.offdiag_pair_j, self.offdiag_pairs_are_exact = _make_offdiag_pairs(
            self.latent_dim,
            int(diagnostic_pair_count),
            int(diagnostic_seed),
        )
        self.offdiag_cross_sum = np.zeros(self.offdiag_pair_i.shape[0], dtype=np.float64)

    def update(self, perturbations: np.ndarray) -> None:
        if perturbations.size == 0:
            return
        block = np.asarray(perturbations, dtype=np.float64)
        self.count += int(block.shape[0])
        self.sum += block.sum(axis=0)
        self.sumsq += np.square(block).sum(axis=0)
        if self.channel_spatial_diagnostics:
            shaped = block.reshape(block.shape[0], *self.latent_shape)
            channel_rows = np.moveaxis(shaped, 1, -1).reshape(-1, self.channel_count)
            self.channel_sample_count += int(channel_rows.shape[0])
            self.channel_sum += channel_rows.sum(axis=0)
            self.channel_outer += channel_rows.T @ channel_rows

            spatial_rows = shaped[:, self.spatial_channel_index, :, :].reshape(block.shape[0], self.spatial_dim)
            self.spatial_sample_count += int(spatial_rows.shape[0])
            self.spatial_sum += spatial_rows.sum(axis=0)
            self.spatial_outer += spatial_rows.T @ spatial_rows
            if self.all_spatial_channel_correlations:
                spatial_by_channel = shaped.reshape(block.shape[0], self.channel_count, self.spatial_dim)
                self.all_spatial_sample_count += int(block.shape[0])
                for channel_index in range(self.channel_count):
                    channel_rows = spatial_by_channel[:, channel_index, :]
                    self.all_spatial_sum[channel_index] += channel_rows.sum(axis=0)
                    self.all_spatial_outer[channel_index] += channel_rows.T @ channel_rows
        local = block[:, : self.local_dims]
        self.local_sum += local.sum(axis=0)
        self.local_outer += local.T @ local
        for offset in range(1, self.diagonal_band_radius + 1):
            width = self.latent_dim - offset
            if width <= 0:
                break
            self.diagonal_band_cross_sum[offset, :width] += (block[:, :width] * block[:, offset:]).sum(axis=0)
        for block_index, start in enumerate(self.diagonal_block_starts):
            stop = int(start) + self.diagonal_block_size
            local_block = block[:, int(start) : stop]
            self.diagonal_block_sum[block_index] += local_block.sum(axis=0)
            self.diagonal_block_outer[block_index] += local_block.T @ local_block
        for start in range(0, self.offdiag_pair_i.size, self.pair_chunk_size):
            stop = min(start + self.pair_chunk_size, self.offdiag_pair_i.size)
            pair_i = self.offdiag_pair_i[start:stop]
            pair_j = self.offdiag_pair_j[start:stop]
            self.offdiag_cross_sum[start:stop] += (block[:, pair_i] * block[:, pair_j]).sum(axis=0)

    def finalize(self) -> dict[str, np.ndarray | int]:
        if self.count < 2:
            raise RuntimeError(f"Need at least two valid NMC pairs, got {self.count}.")
        mean = self.sum / self.count
        variance = (self.sumsq - self.count * np.square(mean)) / (self.count - 1)
        variance = np.nan_to_num(variance, nan=0.0, posinf=0.0, neginf=0.0)
        variance = np.clip(variance, 0.0, None)

        local_mean = self.local_sum / self.count
        local_covariance = (self.local_outer - self.count * np.outer(local_mean, local_mean)) / (self.count - 1)
        local_covariance = np.nan_to_num(local_covariance, nan=0.0, posinf=0.0, neginf=0.0)
        diag = np.sqrt(np.clip(np.diag(local_covariance), 1e-12, None))
        local_correlation = local_covariance / np.outer(diag, diag)
        local_correlation = np.nan_to_num(local_correlation, nan=0.0, posinf=0.0, neginf=0.0)
        eig = np.linalg.eigvalsh(local_covariance.astype(np.float64))
        offdiag_covariance = np.empty((0,), dtype=np.float64)
        offdiag_correlation = np.empty((0,), dtype=np.float64)
        band_offsets = np.arange(-self.diagonal_band_radius, self.diagonal_band_radius + 1, dtype=np.int32)
        band_covariance = np.full((band_offsets.size, self.latent_dim), np.nan, dtype=np.float64)
        band_correlation = np.full_like(band_covariance, np.nan)
        center = self.diagonal_band_radius
        band_covariance[center] = variance
        band_correlation[center] = 1.0
        for offset in range(1, self.diagonal_band_radius + 1):
            width = self.latent_dim - offset
            if width <= 0:
                break
            band_cov = (
                self.diagonal_band_cross_sum[offset, :width]
                - self.count * mean[:width] * mean[offset:]
            ) / (self.count - 1)
            denom = np.sqrt(np.clip(variance[:width] * variance[offset:], 1e-24, None))
            corr_band = np.clip(np.nan_to_num(band_cov / denom, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
            band_covariance[center + offset, :width] = band_cov
            band_correlation[center + offset, :width] = corr_band
            band_covariance[center - offset, offset:] = band_cov
            band_correlation[center - offset, offset:] = corr_band
        block_covariance = np.empty(
            (self.diagonal_block_starts.size, self.diagonal_block_size, self.diagonal_block_size),
            dtype=np.float32,
        )
        block_correlation = np.empty_like(block_covariance)
        for block_index in range(self.diagonal_block_starts.size):
            block_mean = self.diagonal_block_sum[block_index] / self.count
            cov_block = (
                self.diagonal_block_outer[block_index] - self.count * np.outer(block_mean, block_mean)
            ) / (self.count - 1)
            cov_block = np.nan_to_num(cov_block, nan=0.0, posinf=0.0, neginf=0.0)
            diag = np.sqrt(np.clip(np.diag(cov_block), 1e-12, None))
            corr_block = cov_block / np.outer(diag, diag)
            corr_block = np.clip(np.nan_to_num(corr_block, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
            block_covariance[block_index] = cov_block.astype(np.float32)
            block_correlation[block_index] = corr_block.astype(np.float32)
        if self.offdiag_pair_i.size:
            pair_i = self.offdiag_pair_i
            pair_j = self.offdiag_pair_j
            offdiag_covariance = (
                self.offdiag_cross_sum - self.count * mean[pair_i] * mean[pair_j]
            ) / (self.count - 1)
            denom = np.sqrt(np.clip(variance[pair_i] * variance[pair_j], 1e-24, None))
            offdiag_correlation = offdiag_covariance / denom
            offdiag_correlation = np.clip(
                np.nan_to_num(offdiag_correlation, nan=0.0, posinf=0.0, neginf=0.0),
                -1.0,
                1.0,
            )

        channel_correlation = np.empty((0, 0), dtype=np.float64)
        spatial_correlation = np.empty((0, 0), dtype=np.float64)
        spatial_correlations = np.empty((0, 0, 0), dtype=np.float64)
        if self.channel_spatial_diagnostics:
            channel_covariance = _covariance_from_sums(
                self.channel_sum,
                self.channel_outer,
                self.channel_sample_count,
            )
            spatial_covariance = _covariance_from_sums(
                self.spatial_sum,
                self.spatial_outer,
                self.spatial_sample_count,
            )
            channel_correlation = _correlation_from_covariance(channel_covariance)
            spatial_correlation = _correlation_from_covariance(spatial_covariance)
            if self.all_spatial_channel_correlations:
                spatial_correlations = np.empty(
                    (self.channel_count, self.spatial_dim, self.spatial_dim),
                    dtype=np.float32,
                )
                for channel_index in range(self.channel_count):
                    channel_covariance = _covariance_from_sums(
                        self.all_spatial_sum[channel_index],
                        self.all_spatial_outer[channel_index],
                        self.all_spatial_sample_count,
                    )
                    spatial_correlations[channel_index] = _correlation_from_covariance(channel_covariance).astype(np.float32)

        return {
            "perturbation_mean": mean.astype(np.float32),
            "variance": variance.astype(np.float32),
            "B_diag": np.clip(variance, 1e-12, None).astype(np.float32),
            "local_covariance": local_covariance.astype(np.float32),
            "local_correlation": local_correlation.astype(np.float32),
            "local_eigenvalues": eig.astype(np.float32),
            "local_indices": np.arange(self.local_dims, dtype=np.int32),
            "diagonal_band_covariance": band_covariance.astype(np.float32),
            "diagonal_band_correlation": band_correlation.astype(np.float32),
            "diagonal_band_offsets": band_offsets,
            "diagonal_block_covariance": block_covariance,
            "diagonal_block_correlation": block_correlation,
            "diagonal_block_starts": self.diagonal_block_starts,
            "sampled_offdiag_covariance": offdiag_covariance.astype(np.float32),
            "sampled_offdiag_correlation": offdiag_correlation.astype(np.float32),
            "sampled_offdiag_pair_i": self.offdiag_pair_i.astype(np.int32, copy=False),
            "sampled_offdiag_pair_j": self.offdiag_pair_j.astype(np.int32, copy=False),
            "offdiag_pairs_are_exact": bool(self.offdiag_pairs_are_exact),
            "channel_correlation": channel_correlation.astype(np.float32),
            "spatial_correlation": spatial_correlation.astype(np.float32),
            "spatial_correlations": spatial_correlations.astype(np.float32),
            "spatial_channel_index": int(self.spatial_channel_index),
            "valid_pairs": int(self.count),
        }


@torch.inference_mode()
def encode_and_compute_nmc(
    model: CoupledAE,
    dataset,
    device: torch.device,
    batch_size: int,
    max_samples: int,
    lag: int,
    plot_dims: int,
    config: dict[str, Any],
) -> tuple[dict[str, np.ndarray | int], tuple[int, ...], int]:
    if lag < 1:
        raise ValueError("lag must be >= 1")
    loader = build_loader(dataset, batch_size=batch_size, config=config)
    total_samples = min(len(dataset), max_samples)
    latent_shape: tuple[int, ...] | None = None
    stats: OnlineNMCStats | None = None
    history: deque[np.ndarray] = deque(maxlen=lag + 1)
    seen = 0
    start_time = time.perf_counter()
    log_every_batches = int(config.get("log_every_n_batches", 200))
    progress = tqdm(total=total_samples, desc="Encoding NMC latents", unit="sample") if tqdm is not None else None

    LOGGER.info(
        "NMC encoding start: samples=%d lag=%d batch_size=%d loader_workers=%d",
        total_samples,
        lag,
        batch_size,
        int(config.get("num_workers", min(4, os.cpu_count() or 4))),
    )
    for batch_index, batch in enumerate(loader, start=1):
        if seen >= total_samples:
            break
        atm, ocn = _unpack_atm_ocn_batch(batch)
        keep = min(int(atm.shape[0]), total_samples - seen)
        if keep <= 0:
            break
        atm = atm[:keep]
        ocn = ocn[:keep]
        atm = torch.nan_to_num(atm.to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        ocn = torch.nan_to_num(ocn.to(device, non_blocking=True), nan=0.0, posinf=0.0, neginf=0.0)
        with _autocast_context(device, config):
            latent = model.encode(atm, ocn)
        if latent_shape is None:
            latent_shape = tuple(latent.shape[1:])
            latent_dim = int(np.prod(latent_shape))
            stats = OnlineNMCStats(
                latent_dim=latent_dim,
                latent_shape=latent_shape,
                local_dims=plot_dims,
                diagnostic_pair_count=int(config.get("diagnostic_pair_count", 1_000_000)),
                diagnostic_seed=int(config.get("diagnostic_seed", 20260525)),
                pair_chunk_size=int(config.get("diagnostic_pair_chunk_size", 200_000)),
                diagonal_band_radius=int(config.get("diagonal_band_radius", 16)),
                diagonal_block_size=int(config.get("diagonal_block_size", config.get("plot_dims", 1024))),
                diagonal_block_count=int(config.get("diagonal_block_count", 6)),
                diagonal_block_starts=_optional_int_list(config.get("diagonal_block_starts")),
                channel_spatial_diagnostics=bool(config.get("channel_spatial_diagnostics", False)),
                spatial_channel_index=int(config.get("spatial_channel_index", 0)),
                all_spatial_channel_correlations=bool(config.get("plot_all_spatial_channel_correlations", False)),
            )
            LOGGER.info(
                "Latent shape=%s latent_dim=%d local_plot_dims=%d diagonal_band_radius=%d diagonal_block_size=%d diagonal_block_starts=%s diagnostic_pairs=%d exact_pairs=%s channel_spatial_diagnostics=%s spatial_channel=%d all_spatial_channels=%s",
                latent_shape,
                latent_dim,
                stats.local_dims,
                stats.diagonal_band_radius,
                stats.diagonal_block_size,
                stats.diagonal_block_starts.tolist(),
                stats.offdiag_pair_i.size,
                stats.offdiag_pairs_are_exact,
                stats.channel_spatial_diagnostics,
                stats.spatial_channel_index,
                stats.all_spatial_channel_correlations,
            )

        flat = latent.reshape(latent.shape[0], -1).float().cpu().numpy()
        perturb_batch: list[np.ndarray] = []
        for row in flat:
            history.append(row.copy())
            if len(history) > lag:
                perturb_batch.append((history[-1] - history[0]) / np.sqrt(2.0))
        if perturb_batch:
            assert stats is not None
            stats.update(np.stack(perturb_batch, axis=0))

        seen += keep
        if progress is not None:
            progress.update(keep)
        if batch_index % log_every_batches == 0:
            elapsed = max(time.perf_counter() - start_time, 1e-6)
            valid_pairs = 0 if stats is None else stats.count
            LOGGER.info(
                "Encoded %d/%d samples, valid_pairs=%d, rate=%.2f sample/s",
                seen,
                total_samples,
                valid_pairs,
                seen / elapsed,
            )
    if progress is not None:
        progress.close()
    if latent_shape is None or stats is None:
        raise RuntimeError("No samples were encoded for NMC statistics.")
    final = stats.finalize()
    final["sample_count"] = int(seen)
    LOGGER.info(
        "NMC encoding done: samples=%d valid_pairs=%d elapsed=%.1fs",
        seen,
        int(final["valid_pairs"]),
        time.perf_counter() - start_time,
    )
    return final, latent_shape, seen


def plot_bdiag(bdiag: np.ndarray, latent_shape: tuple[int, ...], output_path: Path, cmap: str = "viridis") -> None:
    channels = latent_shape[0]
    field = bdiag.reshape(latent_shape).mean(axis=0) if len(latent_shape) == 3 else bdiag.reshape(1, -1)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    bound = float(np.nanmax(np.abs(field))) if field.size else 1.0
    image = ax.imshow(field, origin="lower", aspect="auto", cmap=cmap, vmin=0.0, vmax=bound)
    ax.set_title(f"NMC B_diag mean over {channels} latent channels")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_matrix(matrix: np.ndarray, output_path: Path, title: str, cmap: str = "bwr") -> None:
    fig, ax = plt.subplots(figsize=(6, 5.5))
    bound = float(np.nanmax(np.abs(matrix))) if matrix.size else 1.0
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, vmin=-bound, vmax=bound)
    ax.set_title(title)
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_correlation_matrix(
    matrix: np.ndarray,
    output_path: Path,
    title: str,
    *,
    cmap: str = "bwr",
    xlabel: str = "Index",
    ylabel: str = "Index",
    dpi: int = PAPER_DPI,
) -> None:
    matrix = np.asarray(matrix, dtype=np.float32)
    if matrix.size == 0:
        LOGGER.warning("Skipping %s because the correlation matrix is empty.", output_path)
        return
    fig_size = (5.6, 4.8) if matrix.shape[0] <= 128 else (6.4, 5.6)
    fig, ax = plt.subplots(figsize=fig_size)
    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0),
        interpolation="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Correlation")
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_paper_correlation_summary(
    channel_correlation: np.ndarray,
    spatial_correlation: np.ndarray,
    output_path: Path,
    *,
    spatial_channel_index: int,
    cmap: str = "bwr",
) -> None:
    channel_correlation = np.asarray(channel_correlation, dtype=np.float32)
    spatial_correlation = np.asarray(spatial_correlation, dtype=np.float32)
    if channel_correlation.size == 0 or spatial_correlation.size == 0:
        LOGGER.warning("Skipping %s because channel/spatial correlations are unavailable.", output_path)
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    axes[0].imshow(channel_correlation, origin="lower", aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    axes[0].set_title("Channel correlation")
    axes[0].set_xlabel("Latent channel")
    axes[0].set_ylabel("Latent channel")

    image = axes[1].imshow(spatial_correlation, origin="lower", aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    axes[1].set_title(f"Spatial correlation (channel {spatial_channel_index})")
    axes[1].set_xlabel("Flattened spatial index")
    axes[1].set_ylabel("Flattened spatial index")
    fig.colorbar(image, ax=axes, fraction=0.035, pad=0.02, label="Correlation")
    fig.savefig(output_path, dpi=PAPER_DPI)
    plt.close(fig)


def _robust_symmetric_limit(values: np.ndarray, percentile: float = 99.5) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 1.0
    limit = float(np.percentile(np.abs(finite), percentile))
    if not np.isfinite(limit) or limit <= 0.0:
        limit = float(np.max(np.abs(finite))) if finite.size else 1.0
    return max(limit, 1e-12)


def plot_local_covariance_block(
    matrix: np.ndarray,
    local_indices: np.ndarray,
    output_path: Path,
    *,
    title: str = "Local covariance block of $B_z$",
    cmap: str = "bwr",
) -> None:
    matrix = np.asarray(matrix)
    start = int(local_indices[0]) if local_indices.size else 0
    end = int(local_indices[-1]) if local_indices.size else matrix.shape[0] - 1
    limit = _robust_symmetric_limit(matrix, percentile=99.7)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)

    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    image = ax.imshow(
        matrix,
        origin="upper",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        rasterized=True,
        extent=[start - 0.5, end + 0.5, end + 0.5, start - 0.5],
    )
    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel("Flattened latent dimension", fontsize=11)
    ax.set_ylabel("Flattened latent dimension", fontsize=11)
    ax.tick_params(axis="both", labelsize=9, length=3)
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label("Covariance", fontsize=10)
    colorbar.ax.tick_params(labelsize=9)
    colorbar.formatter.set_powerlimits((-2, 2))
    colorbar.update_ticks()
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def _make_band_mosaic(band: np.ndarray, tile_width: int, gap: int = 2) -> np.ndarray:
    band = np.asarray(band, dtype=np.float32)
    rows, cols = band.shape
    tile_width = max(int(tile_width), 1)
    n_tiles = int(np.ceil(cols / tile_width))
    mosaic_rows = n_tiles * rows + max(n_tiles - 1, 0) * gap
    mosaic = np.full((mosaic_rows, tile_width), np.nan, dtype=np.float32)
    cursor = 0
    for tile_index in range(n_tiles):
        start = tile_index * tile_width
        stop = min(start + tile_width, cols)
        width = stop - start
        mosaic[cursor : cursor + rows, :width] = band[:, start:stop]
        cursor += rows + gap
    return mosaic


def plot_diagonal_band_mosaic(
    band: np.ndarray,
    offsets: np.ndarray,
    output_path: Path,
    *,
    title: str,
    label: str,
    tile_width: int = 2048,
    cmap: str = "bwr",
    symmetric: bool = True,
) -> None:
    band = np.asarray(band, dtype=np.float32)
    offsets = np.asarray(offsets, dtype=np.int32)
    mosaic = _make_band_mosaic(band, tile_width=tile_width)
    finite = mosaic[np.isfinite(mosaic)]
    if finite.size:
        if symmetric:
            limit = _robust_symmetric_limit(finite, percentile=99.7)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
            vmin = vmax = None
        else:
            norm = None
            vmin = float(np.nanpercentile(finite, 0.5))
            vmax = float(np.nanpercentile(finite, 99.5))
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
                vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
    else:
        norm = None
        vmin, vmax = -1.0, 1.0

    fig_height = float(np.clip(mosaic.shape[0] / 85.0, 5.0, 24.0))
    fig_width = float(np.clip(int(tile_width) / 260.0, 8.0, 18.0))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    image = ax.imshow(
        mosaic,
        origin="upper",
        aspect="auto",
        cmap=cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        rasterized=True,
    )
    rows = band.shape[0]
    n_tiles = int(np.ceil(band.shape[1] / max(int(tile_width), 1)))
    centers = [tile * (rows + 2) + rows / 2 - 0.5 for tile in range(n_tiles)]
    labels = [f"{tile * int(tile_width):,}" for tile in range(n_tiles)]
    max_ticks = 18
    if len(centers) > max_ticks:
        step = int(np.ceil(len(centers) / max_ticks))
        centers = centers[::step]
        labels = labels[::step]
    ax.set_yticks(centers)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"Tile-local latent index, tile width={int(tile_width):,}", fontsize=10)
    ax.set_ylabel("Tile start flattened latent index", fontsize=10)
    ax.set_title(f"{title}\noffsets {int(offsets[0])}..{int(offsets[-1])}", fontsize=12, pad=9)
    ax.tick_params(axis="x", labelsize=8, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.018)
    colorbar.set_label(label, fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    colorbar.formatter.set_powerlimits((-2, 2))
    colorbar.update_ticks()
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_diagonal_block_mosaic(
    blocks: np.ndarray,
    starts: np.ndarray,
    output_path: Path,
    *,
    title: str,
    label: str,
    tile_cols: int = 2,
    cmap: str = "bwr",
) -> None:
    blocks = np.asarray(blocks, dtype=np.float32)
    starts = np.asarray(starts, dtype=np.int32)
    if blocks.size == 0 or blocks.ndim != 3:
        return
    n_blocks, block_size, _ = blocks.shape
    tile_cols = max(1, min(int(tile_cols), n_blocks))
    tile_rows = int(np.ceil(n_blocks / tile_cols))
    gap = max(block_size // 80, 8)
    mosaic = np.full(
        (
            tile_rows * block_size + max(tile_rows - 1, 0) * gap,
            tile_cols * block_size + max(tile_cols - 1, 0) * gap,
        ),
        np.nan,
        dtype=np.float32,
    )
    for block_index in range(n_blocks):
        row = block_index // tile_cols
        col = block_index % tile_cols
        y0 = row * (block_size + gap)
        x0 = col * (block_size + gap)
        mosaic[y0 : y0 + block_size, x0 : x0 + block_size] = blocks[block_index]

    finite = mosaic[np.isfinite(mosaic)]
    limit = _robust_symmetric_limit(finite, percentile=99.7) if finite.size else 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
    fig_width = float(np.clip(mosaic.shape[1] / 360.0, 8.0, 22.0))
    fig_height = float(np.clip(mosaic.shape[0] / 360.0, 6.0, 22.0))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    image = ax.imshow(
        mosaic,
        origin="upper",
        aspect="equal",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        rasterized=True,
    )
    ax.set_title(title, fontsize=12, pad=9)
    ax.set_xlabel("Tile-local flattened latent dimension", fontsize=10)
    ax.set_ylabel("Tile-local flattened latent dimension", fontsize=10)
    ax.tick_params(axis="both", labelsize=8, length=3)
    for block_index, start in enumerate(starts):
        row = block_index // tile_cols
        col = block_index % tile_cols
        y0 = row * (block_size + gap)
        x0 = col * (block_size + gap)
        ax.text(
            x0 + 0.02 * block_size,
            y0 + 0.06 * block_size,
            f"start={int(start):,}",
            color="#222222",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.78},
        )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.018)
    colorbar.set_label(label, fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    colorbar.formatter.set_powerlimits((-2, 2))
    colorbar.update_ticks()
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_diagonal_block_tiles(
    blocks: np.ndarray,
    starts: np.ndarray,
    output_dir: Path,
    *,
    prefix: str,
    title: str,
    label: str,
    cmap: str = "bwr",
) -> None:
    blocks = np.asarray(blocks, dtype=np.float32)
    starts = np.asarray(starts, dtype=np.int32)
    if blocks.size == 0 or blocks.ndim != 3:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    finite = blocks[np.isfinite(blocks)]
    limit = _robust_symmetric_limit(finite, percentile=99.7) if finite.size else 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
    for block_index, start in enumerate(starts):
        block = blocks[block_index]
        fig, ax = plt.subplots(figsize=(7.2, 6.3), constrained_layout=True)
        image = ax.imshow(
            block,
            origin="upper",
            aspect="equal",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )
        end = int(start) + block.shape[0] - 1
        ax.set_title(f"{title}\nflattened latent {int(start):,}:{end:,}", fontsize=12, pad=9)
        ax.set_xlabel("Flattened latent dimension", fontsize=10)
        ax.set_ylabel("Flattened latent dimension", fontsize=10)
        ax.tick_params(axis="both", labelsize=8, length=3)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
        colorbar.set_label(label, fontsize=10)
        colorbar.ax.tick_params(labelsize=8)
        colorbar.formatter.set_powerlimits((-2, 2))
        colorbar.update_ticks()
        fig.savefig(output_dir / f"{prefix}_{block_index:03d}_start_{int(start):06d}.png", dpi=PAPER_DPI, bbox_inches="tight")
        plt.close(fig)


def plot_variance_histogram(variance: np.ndarray, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(np.asarray(variance).reshape(-1), bins=40, color="#2f6c8f", edgecolor="white")
    ax.set_title("Histogram of Latent Background Variance")
    ax.set_xlabel("Variance")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def percentile_from_histogram(hist: np.ndarray, bin_edges: np.ndarray, percentile: float) -> float:
    total = int(hist.sum())
    if total <= 0:
        return 0.0
    target = percentile / 100.0 * total
    cumulative = np.cumsum(hist)
    idx = int(np.searchsorted(cumulative, target, side="left"))
    idx = min(max(idx, 0), hist.size - 1)
    prev = float(cumulative[idx - 1]) if idx > 0 else 0.0
    count = float(hist[idx])
    if count <= 0:
        return float(0.5 * (bin_edges[idx] + bin_edges[idx + 1]))
    fraction = np.clip((target - prev) / count, 0.0, 1.0)
    return float(bin_edges[idx] + fraction * (bin_edges[idx + 1] - bin_edges[idx]))


def compute_offdiag_correlation_distribution(
    offdiag_correlation: np.ndarray,
    corr_bins: int,
    mode: str,
) -> dict[str, np.ndarray | float | int | str]:
    corr_bins = max(3, int(corr_bins))
    signed_edges = np.linspace(-1.0, 1.0, corr_bins + 1, dtype=np.float64)
    abs_edges = np.linspace(0.0, 1.0, corr_bins + 1, dtype=np.float64)
    values = np.asarray(offdiag_correlation, dtype=np.float64).reshape(-1)
    values = values[np.isfinite(values)]
    values = np.clip(values, -1.0, 1.0)
    signed_hist = np.histogram(values, bins=signed_edges)[0].astype(np.int64)
    abs_hist = np.histogram(np.abs(values), bins=abs_edges)[0].astype(np.int64)
    return {
        "offdiag_corr_hist": signed_hist,
        "offdiag_corr_bin_edges": signed_edges,
        "offdiag_abs_corr_hist": abs_hist,
        "offdiag_abs_corr_bin_edges": abs_edges,
        "offdiag_corr_p95": percentile_from_histogram(signed_hist, signed_edges, 95.0),
        "offdiag_abs_corr_p95": percentile_from_histogram(abs_hist, abs_edges, 95.0),
        "offdiag_corr_pair_count": int(values.size),
        "offdiag_corr_distribution_mode": mode,
    }


def compute_sampled_diag_approx_diagnostics(
    variance: np.ndarray,
    sampled_offdiag_covariance: np.ndarray,
    sampled_offdiag_correlation: np.ndarray,
    latent_dim: int,
) -> dict[str, float]:
    variance = np.asarray(variance, dtype=np.float64).reshape(-1)
    offdiag_covariance = np.asarray(sampled_offdiag_covariance, dtype=np.float64).reshape(-1)
    offdiag_correlation = np.asarray(sampled_offdiag_correlation, dtype=np.float64).reshape(-1)
    offdiag_correlation = offdiag_correlation[np.isfinite(offdiag_correlation)]
    offdiag_abs_corr = np.abs(offdiag_correlation)
    diag_energy = float(np.square(variance).sum())
    if offdiag_covariance.size:
        mean_offdiag_sq = float(np.nanmean(np.square(offdiag_covariance)))
        offdiag_energy_estimate = mean_offdiag_sq * float(latent_dim) * float(max(latent_dim - 1, 0))
    else:
        offdiag_energy_estimate = 0.0
    total = diag_energy + offdiag_energy_estimate
    return {
        "diag_energy": diag_energy,
        "offdiag_energy_estimate": offdiag_energy_estimate,
        "diag_energy_fraction_estimate": float(diag_energy / total) if total > 0 else 0.0,
        "offdiag_energy_fraction_estimate": float(offdiag_energy_estimate / total) if total > 0 else 0.0,
        "offdiag_abs_corr_median": float(np.median(offdiag_abs_corr)) if offdiag_abs_corr.size else 0.0,
        "offdiag_abs_corr_p95": float(np.percentile(offdiag_abs_corr, 95)) if offdiag_abs_corr.size else 0.0,
    }


def plot_offdiag_correlation_distribution(
    corr_distribution: dict[str, np.ndarray | float | int | str],
    output_path: Path,
    *,
    color: str = "#1f77b4",
    p95_color: str = "#d88700",
) -> None:
    hist = np.asarray(corr_distribution["offdiag_corr_hist"], dtype=np.float64)
    bin_edges = np.asarray(corr_distribution["offdiag_corr_bin_edges"], dtype=np.float64)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    density = hist / max(float(hist.sum()), 1.0) / widths
    abs_p95 = float(corr_distribution["offdiag_abs_corr_p95"])

    fig, ax = plt.subplots(figsize=(7.8, 4.9), constrained_layout=True)
    ax.plot(centers, density, color=color, linewidth=2.2, solid_capstyle="round")
    ax.axvline(0.0, color="#4c4c4c", linewidth=1.1, alpha=0.75)
    ax.axvline(abs_p95, color=p95_color, linewidth=1.9, linestyle="--", label=f"95% |corr| = {abs_p95:.3f}")
    ax.axvline(-abs_p95, color=p95_color, linewidth=1.9, linestyle="--")
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_title("Off-diagonal correlation distribution", fontsize=13, pad=10)
    ax.set_xlabel(r"$\mathrm{corr}(i,j),\ i \ne j$", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.tick_params(axis="both", labelsize=9, length=3)
    ax.grid(True, color="#c7c7c7", linewidth=0.8, alpha=0.35)
    ax.legend()
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.savefig(output_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def compute_diag_approx_diagnostics(
    local_covariance: np.ndarray,
    local_correlation: np.ndarray,
) -> dict[str, np.ndarray | float]:
    if local_covariance.shape[0] != local_covariance.shape[1]:
        raise ValueError("local_covariance must be square.")
    offdiag_mask = ~np.eye(local_covariance.shape[0], dtype=bool)
    offdiag_abs_corr = np.abs(local_correlation[offdiag_mask])
    diag_entries = np.diag(local_covariance).astype(np.float64, copy=False)
    offdiag_entries = local_covariance.astype(np.float64, copy=True)
    np.fill_diagonal(offdiag_entries, 0.0)
    diag_energy = float(np.square(diag_entries).sum())
    offdiag_energy = float(np.square(offdiag_entries).sum())
    total_energy = diag_energy + offdiag_energy
    row_coupling_ratio = np.abs(offdiag_entries).sum(axis=1) / (np.abs(diag_entries) + 1e-12)
    return {
        "offdiag_abs_corr": offdiag_abs_corr.astype(np.float32),
        "diag_energy": diag_energy,
        "offdiag_energy": offdiag_energy,
        "diag_energy_fraction": float(diag_energy / total_energy) if total_energy > 0 else 0.0,
        "offdiag_energy_fraction": float(offdiag_energy / total_energy) if total_energy > 0 else 0.0,
        "row_coupling_ratio": row_coupling_ratio.astype(np.float32),
        "offdiag_abs_corr_median": float(np.median(offdiag_abs_corr)) if offdiag_abs_corr.size else 0.0,
        "offdiag_abs_corr_p95": float(np.percentile(offdiag_abs_corr, 95)) if offdiag_abs_corr.size else 0.0,
        "row_coupling_ratio_median": float(np.median(row_coupling_ratio)) if row_coupling_ratio.size else 0.0,
        "row_coupling_ratio_p95": float(np.percentile(row_coupling_ratio, 95)) if row_coupling_ratio.size else 0.0,
    }


def plot_offdiag_corr_histogram(offdiag_abs_corr: np.ndarray, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(offdiag_abs_corr, bins=50, color="#3f7f93", edgecolor="white")
    if offdiag_abs_corr.size:
        median = float(np.median(offdiag_abs_corr))
        p95 = float(np.percentile(offdiag_abs_corr, 95))
        ax.axvline(median, color="#a23b72", linewidth=2, label=f"median={median:.3g}")
        ax.axvline(p95, color="#f18f01", linewidth=2, linestyle="--", label=f"p95={p95:.3g}")
        ax.legend()
    ax.set_title("Absolute Off-Diagonal Correlation Histogram")
    ax.set_xlabel("|corr(i, j)|, i != j")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_offdiag_corr_cdf(offdiag_abs_corr: np.ndarray, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if offdiag_abs_corr.size:
        sorted_vals = np.sort(offdiag_abs_corr)
        cdf = np.linspace(0.0, 1.0, sorted_vals.size, endpoint=True)
        ax.plot(sorted_vals, cdf, linewidth=2.0, color="#1b6ca8")
    ax.set_title("CDF of Absolute Off-Diagonal Correlation")
    ax.set_xlabel("|corr(i, j)|, i != j")
    ax.set_ylabel("CDF")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_diag_energy_fraction(diag_stats: dict[str, np.ndarray | float], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    diag_fraction = float(diag_stats["diag_energy_fraction"])
    offdiag_fraction = float(diag_stats["offdiag_energy_fraction"])
    ax.bar(["diag", "offdiag"], [diag_fraction, offdiag_fraction], color=["#2a9d8f", "#e76f51"])
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Frobenius Energy Fraction in Local Block")
    ax.set_ylabel("Energy fraction")
    for index, value in enumerate([diag_fraction, offdiag_fraction]):
        ax.text(index, min(value + 0.02, 0.98), f"{value:.3f}", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_row_coupling_ratio_histogram(diag_stats: dict[str, np.ndarray | float], output_path: Path) -> None:
    row_coupling_ratio = np.asarray(diag_stats["row_coupling_ratio"])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(row_coupling_ratio, bins=40, color="#577590", edgecolor="white")
    if row_coupling_ratio.size:
        median = float(diag_stats["row_coupling_ratio_median"])
        p95 = float(diag_stats["row_coupling_ratio_p95"])
        ax.axvline(median, color="#a23b72", linewidth=2, label=f"median={median:.3g}")
        ax.axvline(p95, color="#f18f01", linewidth=2, linestyle="--", label=f"p95={p95:.3g}")
        ax.legend()
    ax.set_title("Row-Wise Offdiag/Diag Strength Ratio")
    ax.set_xlabel("sum_j!=i |B_ij| / |B_ii|")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run(config: dict) -> Path:
    device = resolve_device(config.get("device"))
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = Path(config["data_path"])
    max_samples = config.get("max_samples")
    if max_samples is None:
        max_samples = DEFAULT_MAX_SAMPLES
    max_samples = int(max_samples)
    split = str(config.get("split", config.get("dataset_split", "val"))).lower()
    val_fraction = float(config.get("val_fraction", 0.1))
    if data_path.is_dir() and (data_path / FRAME_MANIFEST).exists():
        dataset = CM2FrameDataset(
            data_path,
            split=split,
            val_fraction=val_fraction,
            max_samples=max_samples,
            load_mmap=bool(config.get("mmap_frames", False)),
        )
    else:
        dataset = CM2TensorDataset(data_path, metadata_path=config.get("metadata_path"), normalize=True)
        if split != "all":
            total = len(dataset)
            n_val = max(1, int(total * val_fraction)) if total > 1 else 0
            order = list(range(total))
            indices = order[-n_val:] if split == "val" and n_val > 0 else order[:-n_val] if n_val > 0 else order
            dataset = CM2DatasetView(dataset, indices[:max_samples])
    LOGGER.info(
        "Dataset ready: path=%s split=%s val_fraction=%.3f samples=%d max_samples=%d device=%s output=%s",
        data_path,
        split,
        val_fraction,
        len(dataset),
        max_samples,
        device,
        output_dir,
    )
    model = load_model(config["checkpoint"], device)
    LOGGER.info("Checkpoint loaded: %s", config["checkpoint"])
    stats, latent_shape, encoded_samples = encode_and_compute_nmc(
        model,
        dataset,
        device=device,
        batch_size=int(config.get("batch_size", 1)),
        max_samples=max_samples,
        lag=int(config.get("lag", 1)),
        plot_dims=int(config.get("plot_dims", 512)),
        config=config,
    )
    time_values = list(getattr(dataset, "time_values", list(range(len(dataset)))))[:encoded_samples]
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
    LOGGER.info("Saving NMC arrays to %s", npz_path)
    np.savez_compressed(
        npz_path,
        variance=stats["variance"],
        B_diag=stats["B_diag"],
        background_covariance_diag=stats["B_diag"],
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
        lag=np.array([int(config.get("lag", 1))], dtype=np.int32),
        split=np.array([split], dtype=object),
        val_fraction=np.array([val_fraction], dtype=np.float64),
        diagnostic_pair_count=np.array([corr_distribution["offdiag_corr_pair_count"]], dtype=np.int64),
        diagnostic_seed=np.array([int(config.get("diagnostic_seed", 20260525))], dtype=np.int64),
        diagnostic_pairs_are_exact=np.array([bool(stats.get("offdiag_pairs_are_exact", False))], dtype=bool),
        diagonal_band_radius=np.array([int(config.get("diagonal_band_radius", 16))], dtype=np.int32),
        diagonal_block_size=np.array([int(config.get("diagonal_block_size", config.get("plot_dims", 1024)))], dtype=np.int32),
        diagonal_block_starts=stats["diagonal_block_starts"],
        channel_correlation=stats["channel_correlation"],
        spatial_correlation=stats["spatial_correlation"],
        spatial_correlations=stats["spatial_correlations"],
        spatial_channel_index=np.array([stats["spatial_channel_index"]], dtype=np.int32),
        sample_count=np.array([stats["sample_count"]], dtype=np.int32),
        valid_pairs=np.array([stats["valid_pairs"]], dtype=np.int32),
        time_values=np.array(time_values, dtype=object),
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
    LOGGER.info("Writing NMC figures")
    plot_bdiag(
        stats["B_diag"],
        latent_shape,
        output_dir / "latent_nmc_Bz_diag.png",
        cmap=str(config.get("diag_cmap", "viridis")),
    )
    if bool(config.get("plot_local_block", True)):
        plot_local_covariance_block(
            stats["local_covariance"],
            stats["local_indices"],
            output_dir / "latent_nmc_Bz_local_block.png",
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
    if bool(config.get("plot_diagonal_band", True)):
        diagonal_tile_width = int(config.get("diagonal_band_tile_width", 1024))
        plot_diagonal_band_mosaic(
            stats["diagonal_band_covariance"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_Bz_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of NMC $B_z$",
            label="Covariance",
            tile_width=diagonal_tile_width,
            cmap=str(config.get("matrix_cmap", "bwr")),
            symmetric=True,
        )
        plot_diagonal_band_mosaic(
            stats["diagonal_band_correlation"],
            stats["diagonal_band_offsets"],
            output_dir / "latent_nmc_corr_diagonal_band_mosaic.png",
            title="Full diagonal band mosaic of NMC correlation",
            label="Correlation",
            tile_width=diagonal_tile_width,
            cmap=str(config.get("matrix_cmap", "bwr")),
            symmetric=True,
        )
    if bool(config.get("plot_diagonal_blocks", True)) and stats["diagonal_block_covariance"].size:
        block_tile_cols = int(config.get("diagonal_block_tile_cols", 2))
        plot_diagonal_block_mosaic(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            output_dir / "latent_nmc_Bz_diagonal_block_mosaic.png",
            title=f"Diagonal {int(config.get('diagonal_block_size', config.get('plot_dims', 1024))):,}x{int(config.get('diagonal_block_size', config.get('plot_dims', 1024))):,} blocks of NMC $B_z$",
            label="Covariance",
            tile_cols=block_tile_cols,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_mosaic(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            output_dir / "latent_nmc_corr_diagonal_block_mosaic.png",
            title=f"Diagonal {int(config.get('diagonal_block_size', config.get('plot_dims', 1024))):,}x{int(config.get('diagonal_block_size', config.get('plot_dims', 1024))):,} blocks of NMC correlation",
            label="Correlation",
            tile_cols=block_tile_cols,
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    if bool(config.get("save_diagonal_block_tiles", True)) and stats["diagonal_block_covariance"].size:
        block_dir = output_dir / "diagonal_block_tiles"
        plot_diagonal_block_tiles(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="Bz_covariance",
            title="NMC $B_z$ diagonal block",
            label="Covariance",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
        plot_diagonal_block_tiles(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="correlation",
            title="NMC correlation diagonal block",
            label="Correlation",
            cmap=str(config.get("matrix_cmap", "bwr")),
        )
    if bool(config.get("plot_local_block", True)):
        plot_matrix(
            stats["local_correlation"],
            output_dir / "latent_nmc_local_correlation.png",
            "Local NMC Correlation Block",
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
    LOGGER.info("NMC output complete: %s", npz_path)
    return npz_path


def _safe_name(value: str) -> str:
    name = Path(value).stem
    name = re.sub(r"^checkpoint_", "ckpt_", name)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return name or "nmc_version"


def _checkpoint_entries(config: dict, cli_checkpoints: list[str] | None, cli_labels: list[str] | None) -> list[tuple[str, str | None]]:
    if cli_checkpoints:
        labels = cli_labels or []
        return [
            (checkpoint, labels[index] if index < len(labels) else None)
            for index, checkpoint in enumerate(cli_checkpoints)
        ]

    config_checkpoints = config.get("checkpoints")
    if isinstance(config_checkpoints, list) and config_checkpoints:
        entries: list[tuple[str, str | None]] = []
        for item in config_checkpoints:
            if isinstance(item, dict):
                path = item.get("path") or item.get("checkpoint")
                if path:
                    entries.append((str(path), item.get("label") or item.get("name")))
            else:
                entries.append((str(item), None))
        if entries:
            return entries

    return [(str(config["checkpoint"]), None)]


def run_versions(config: dict, checkpoints: list[str] | None, labels: list[str] | None) -> list[Path]:
    entries = _checkpoint_entries(config, checkpoints, labels)
    if len(entries) == 1:
        run_config = dict(config)
        run_config["checkpoint"] = entries[0][0]
        return [run(run_config)]

    base_output_dir = Path(config["output_dir"])
    outputs: list[Path] = []
    for checkpoint, label in entries:
        run_config = dict(config)
        version_name = label or _safe_name(checkpoint)
        run_config["checkpoint"] = checkpoint
        run_config["output_dir"] = str(base_output_dir / version_name)
        LOGGER.info("Running NMC version %s -> %s", version_name, run_config["output_dir"])
        outputs.append(run(run_config))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "nmc_background.json"))
    parser.add_argument("--checkpoint", nargs="+", default=None, help="Override config checkpoint. Pass multiple paths to write one subdir per version.")
    parser.add_argument("--labels", nargs="+", default=None, help="Optional labels for multiple --checkpoint paths.")
    parser.add_argument("--output-dir", default=None, help="Override config output_dir.")
    args = parser.parse_args()
    configure_logging()
    with open(args.config, "r", encoding="utf-8") as f:
        config = json.load(f)
    if args.output_dir:
        config["output_dir"] = args.output_dir
    outputs = run_versions(config, args.checkpoint, args.labels)
    print("saved " + ", ".join(str(path) for path in outputs), flush=True)


if __name__ == "__main__":
    main()
