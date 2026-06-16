from __future__ import annotations

import argparse
from collections import Counter, deque
from datetime import timedelta
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import torch
from torch.utils.data import DataLoader, Subset
from tqdm import tqdm


PAPER_DPI = 300

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from dataset import WRFSlidingDataset
from networks.transformer import LGUnet_all


DATA_CONFIG = {
    "nc_dir": "/data/zsq_data/d01_d02_d03/d03_wrf",
    "input_len": 1,
    "label_len": 1,
    "input_vars": ["U10", "V10", "T2", "Q2"],
    "label_vars": ["U10", "V10", "T2", "Q2"],
}

MODEL_CONFIG = {
    "img_size": [256, 256],
    "patch_size": 4,
    "stride": [4, 4],
    "in_chans": 4,
    "out_chans": 4,
    "enc_depths": [2, 2, 6],
    "enc_heads": [3, 6, 12],
    "lg_depths": [2, 2],
    "lg_heads": [6, 12],
    "inchans_list": [4],
    "outchans_list": [4],
    "enc_dim": 96,
    "embed_dim": 768,
    "window_size": 8,
    "Weather_T": 1,
    "use_checkpoint": False,
    "pre_norm": True,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use quarter 11+12 samples to estimate latent-space background covariance with NMC."
    )
    parser.add_argument("--nc-dir", default=DATA_CONFIG["nc_dir"], help="NetCDF directory.")
    parser.add_argument("--model-path", default=None, help="Checkpoint path. Defaults to best_model.pth if it exists.")
    parser.add_argument("--output-dir", default=None, help="Directory for npz and figures.")
    parser.add_argument("--base-year", type=int, default=2022, help="Base year used by get_sample_quarters().")
    parser.add_argument("--quarters", type=int, nargs="+", default=[11, 12], help="Quarter ids used for the statistics set. Defaults to val+test quarters 11 and 12.")
    parser.add_argument("--lag", type=int, default=1, help="Temporal lag for NMC differences in sample steps.")
    parser.add_argument("--batch-size", type=int, default=2, help="Batch size for latent extraction.")
    parser.add_argument("--num-workers", type=int, default=0, help="DataLoader workers.")
    parser.add_argument("--max-samples", type=int, default=5000, help="Optional cap for the statistics-set size.")
    parser.add_argument("--plot-dims", type=int, default=1, help="Number of flattened latent dimensions shown in the optional local covariance plot.")
    parser.add_argument("--plot-offset", type=int, default=0, help="Start index of the flattened latent block shown in the local covariance plot.")
    parser.add_argument("--plot-local-block", action=argparse.BooleanOptionalAction, default=False, help="Plot the local flattened covariance block.")
    parser.add_argument("--corr-bins", type=int, default=201, help="Bins used for the off-diagonal correlation distribution.")
    parser.add_argument("--diagnostic-pair-count", type=int, default=1_000_000, help="Global non-diagonal latent pairs used for off-diagonal diagnostics.")
    parser.add_argument("--diagnostic-pair-batch-size", type=int, default=200_000, help="Pair batch size for sampled off-diagonal diagnostics.")
    parser.add_argument("--diagnostic-seed", type=int, default=20260525, help="Random seed for sampled off-diagonal latent pairs.")
    parser.add_argument("--channel-spatial-diagnostics", action=argparse.BooleanOptionalAction, default=True, help="Estimate channel-level and spatial latent correlation matrices.")
    parser.add_argument("--plot-channel-spatial-correlation", action=argparse.BooleanOptionalAction, default=True, help="Plot channel-level and spatial latent correlation matrices.")
    parser.add_argument("--plot-all-spatial-channel-correlations", action=argparse.BooleanOptionalAction, default=True, help="Plot one spatial correlation matrix for every latent channel.")
    parser.add_argument("--latent-channel-axis", type=int, default=-1, help="Channel axis inside latent_shape. WRF latents are usually H,W,C, so the default is -1.")
    parser.add_argument("--spatial-channel-index", type=int, default=0, help="Latent channel used for the spatial correlation matrix.")
    parser.add_argument("--diagonal-band-radius", type=int, default=16, help="Half-width around the full flattened latent diagonal for narrow-band mosaic plots.")
    parser.add_argument("--diagonal-band-tile-width", type=int, default=1024, help="Tile width for full diagonal narrow-band mosaic plots.")
    parser.add_argument("--plot-diagonal-band", action=argparse.BooleanOptionalAction, default=False, help="Compute and plot the flattened diagonal narrow-band mosaic.")
    parser.add_argument("--diagonal-block-size", type=int, default=2048, help="Size of each diagonal covariance/correlation block.")
    parser.add_argument("--diagonal-block-count", type=int, default=0, help="Number of diagonal blocks sampled along the full flattened latent diagonal.")
    parser.add_argument("--diagonal-block-starts", type=int, nargs="+", default=None, help="Optional explicit flattened latent starts for diagonal block plots.")
    parser.add_argument("--diagonal-block-tile-cols", type=int, default=4, help="Number of columns in diagonal block mosaic figures.")
    parser.add_argument("--plot-diagonal-blocks", action=argparse.BooleanOptionalAction, default=False, help="Compute and plot large flattened diagonal covariance/correlation blocks.")
    parser.add_argument("--save-diagonal-block-tiles", action=argparse.BooleanOptionalAction, default=False, help="Save every diagonal block as its own figure.")
    parser.add_argument("--corr-sample-pairs", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--corr-pair-batch-size", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--corr-random-seed", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--device", default=None, help="cuda / cpu. Defaults to cuda when available.")
    return parser.parse_args()


def build_model() -> LGUnet_all:
    return LGUnet_all(**MODEL_CONFIG)


def resolve_model_path(requested_path: str | None) -> Path:
    if requested_path:
        model_path = Path(requested_path).expanduser().resolve()
        if not model_path.exists():
            raise FileNotFoundError(f"Model checkpoint not found: {model_path}")
        return model_path

    base_dir = Path(__file__).resolve().parent
    candidates = [
        base_dir / "best_model.pth",
        base_dir / "model_forgoal.pth",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError("No checkpoint found. Please provide --model-path.")


def load_checkpoint(model: torch.nn.Module, model_path: Path, device: torch.device) -> None:
    checkpoint = torch.load(model_path, map_location=device)
    if any(key.startswith("module.") for key in checkpoint.keys()):
        checkpoint = {key[7:]: value for key, value in checkpoint.items()}
    model.load_state_dict(checkpoint)


def resolve_quarters(dataset: WRFSlidingDataset, base_year: int, quarters: list[int] | None) -> tuple[list[int], list[int]]:
    quarter_ids = dataset.get_sample_quarters(base_year=base_year)
    available_quarters = sorted(set(quarter_ids))
    if quarters is None or len(quarters) == 0:
        return available_quarters, quarter_ids
    return quarters, quarter_ids


def build_stats_subset(dataset: WRFSlidingDataset, base_year: int, quarters: list[int] | None, max_samples: int | None) -> tuple[Subset, list[int], list, list[int]]:
    selected_quarters, quarter_ids = resolve_quarters(dataset, base_year, quarters)
    stats_indices = [idx for idx, quarter in enumerate(quarter_ids) if quarter in selected_quarters]
    if not stats_indices:
        raise RuntimeError(f"No samples found for quarters {selected_quarters}.")

    if max_samples is not None:
        stats_indices = stats_indices[:max_samples]

    subset = Subset(dataset, stats_indices)
    subset_times = [dataset.sample_start_times[idx] for idx in stats_indices]
    return subset, stats_indices, subset_times, selected_quarters


def infer_sample_step(sample_times: list) -> timedelta:
    if len(sample_times) < 2:
        raise RuntimeError("At least two samples are required to infer the temporal step.")

    deltas = [sample_times[idx + 1] - sample_times[idx] for idx in range(len(sample_times) - 1)]
    positive_deltas = [delta for delta in deltas if delta.total_seconds() > 0]
    if not positive_deltas:
        raise RuntimeError("Failed to infer a positive temporal step from the selected samples.")

    return Counter(positive_deltas).most_common(1)[0][0]


def flatten_latent_batch(latent_batch: torch.Tensor) -> tuple[torch.Tensor, tuple[int, ...]]:
    latent_shape = tuple(latent_batch.shape[1:])
    flat_latent = latent_batch.reshape(latent_batch.shape[0], -1)
    return flat_latent, latent_shape


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


def _correlation_from_covariance(covariance: np.ndarray) -> np.ndarray:
    covariance = np.asarray(covariance, dtype=np.float64)
    if covariance.size == 0:
        return covariance.astype(np.float32)
    std = np.sqrt(np.clip(np.diag(covariance), 1e-24, None))
    correlation = covariance / np.outer(std, std)
    correlation = np.clip(np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
    return correlation.astype(np.float32)


def _compute_channel_spatial_correlations(
    perturbation_centered: np.ndarray,
    latent_shape: tuple[int, ...],
    channel_axis: int,
    spatial_channel_index: int,
    all_spatial_channel_correlations: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    if len(latent_shape) < 2:
        return (
            np.empty((0, 0), dtype=np.float32),
            np.empty((0, 0), dtype=np.float32),
            np.empty((0, 0, 0), dtype=np.float32),
            0,
            0,
        )
    if int(np.prod(latent_shape)) != int(perturbation_centered.shape[1]):
        raise ValueError("latent_shape does not match flattened perturbation dimension.")

    normalized_channel_axis = int(channel_axis) % len(latent_shape)
    channel_count = int(latent_shape[normalized_channel_axis])
    selected_channel = min(max(int(spatial_channel_index), 0), channel_count - 1)

    shaped = perturbation_centered.reshape((perturbation_centered.shape[0],) + tuple(latent_shape))
    channel_last = np.moveaxis(shaped, 1 + normalized_channel_axis, -1)
    channel_rows = channel_last.reshape(-1, channel_count)
    channel_covariance = (channel_rows.T @ channel_rows) / max(channel_rows.shape[0] - 1, 1)
    channel_correlation = _correlation_from_covariance(channel_covariance)

    spatial_values = np.take(shaped, selected_channel, axis=1 + normalized_channel_axis)
    spatial_rows = spatial_values.reshape(perturbation_centered.shape[0], -1)
    spatial_covariance = (spatial_rows.T @ spatial_rows) / max(spatial_rows.shape[0] - 1, 1)
    spatial_correlation = _correlation_from_covariance(spatial_covariance)

    if all_spatial_channel_correlations:
        spatial_dim = int(spatial_rows.shape[1])
        spatial_correlations = np.empty((channel_count, spatial_dim, spatial_dim), dtype=np.float32)
        for channel_index in tqdm(range(channel_count), desc="Spatial channel correlations", unit="channel"):
            channel_spatial_rows = np.take(shaped, channel_index, axis=1 + normalized_channel_axis).reshape(
                perturbation_centered.shape[0],
                spatial_dim,
            )
            channel_spatial_covariance = (channel_spatial_rows.T @ channel_spatial_rows) / max(channel_spatial_rows.shape[0] - 1, 1)
            spatial_correlations[channel_index] = _correlation_from_covariance(channel_spatial_covariance)
    else:
        spatial_correlations = np.empty((0, 0, 0), dtype=np.float32)

    return channel_correlation, spatial_correlation, spatial_correlations, normalized_channel_axis, selected_channel


@torch.inference_mode()
def collect_nmc_statistics(
    model: torch.nn.Module,
    dataloader: DataLoader,
    sample_times: list,
    device: torch.device,
    lag: int,
    plot_dims: int,
    plot_offset: int,
    diagonal_band_radius: int,
    diagonal_block_size: int,
    diagonal_block_count: int,
    diagonal_block_starts: list[int] | None,
    compute_diagonal_band: bool,
    compute_diagonal_blocks: bool,
    channel_spatial_diagnostics: bool,
    latent_channel_axis: int,
    spatial_channel_index: int,
    all_spatial_channel_correlations: bool,
) -> dict[str, np.ndarray | int | float]:
    expected_step = infer_sample_step(sample_times)
    expected_pair_step = expected_step * lag

    latent_batches: list[np.ndarray] = []
    latent_shape: tuple[int, ...] | None = None

    for inputs, _, _ in tqdm(dataloader, desc="Encoding latent set"):
        inputs = inputs.to(device, non_blocking=(device.type == "cuda"))
        latent_batch = model.encode(inputs)
        flat_latent, current_shape = flatten_latent_batch(latent_batch)
        if latent_shape is None:
            latent_shape = current_shape
        latent_batches.append(flat_latent.detach().cpu().numpy().astype(np.float32))

    if not latent_batches:
        raise RuntimeError("No latent vectors were extracted.")

    latent_matrix = np.concatenate(latent_batches, axis=0)
    if latent_shape is None or latent_matrix.shape[0] < 2:
        raise RuntimeError("The statistics set is too small to estimate covariance.")

    perturbations = []
    valid_pairs = 0
    for idx in range(lag, latent_matrix.shape[0]):
        if sample_times[idx] - sample_times[idx - lag] == expected_pair_step:
            perturbations.append((latent_matrix[idx] - latent_matrix[idx - lag]) / np.sqrt(2.0))
            valid_pairs += 1

    if valid_pairs < 2:
        raise RuntimeError("Valid NMC pairs are fewer than 2. Try a smaller --lag or check whether the time series is continuous.")

    perturbation_matrix = np.stack(perturbations, axis=0).astype(np.float64, copy=False)
    perturbation_mean = perturbation_matrix.mean(axis=0)
    perturbation_centered = perturbation_matrix - perturbation_mean

    feature_dim = perturbation_matrix.shape[1]
    local_start = max(0, plot_offset)
    local_end = min(feature_dim, local_start + max(1, plot_dims))
    if local_start >= local_end:
        local_start = 0
        local_end = min(feature_dim, max(1, plot_dims))
    local_indices = np.arange(local_start, local_end, dtype=np.int32)

    variance = np.square(perturbation_centered.astype(np.float64, copy=False)).sum(axis=0, dtype=np.float64) / (valid_pairs - 1)
    local_centered = perturbation_centered[:, local_start:local_end]
    local_covariance = (local_centered.T @ local_centered) / (valid_pairs - 1)
    local_std = np.sqrt(np.clip(np.diag(local_covariance), a_min=1e-12, a_max=None))
    local_correlation = local_covariance / np.outer(local_std, local_std)
    local_correlation = np.nan_to_num(local_correlation, nan=0.0, posinf=0.0, neginf=0.0)
    local_eigenvalues = np.linalg.eigvalsh(local_covariance)[::-1]
    if compute_diagonal_band:
        diagonal_band_radius = max(int(diagonal_band_radius), 0)
        band_offsets = np.arange(-diagonal_band_radius, diagonal_band_radius + 1, dtype=np.int32)
        diagonal_band_covariance = np.full((band_offsets.size, feature_dim), np.nan, dtype=np.float32)
        diagonal_band_correlation = np.full_like(diagonal_band_covariance, np.nan)
        center = diagonal_band_radius
        diagonal_band_covariance[center] = variance.astype(np.float32)
        diagonal_band_correlation[center] = 1.0
        for offset in range(1, diagonal_band_radius + 1):
            width = feature_dim - offset
            if width <= 0:
                break
            cov_band = np.einsum(
                "ni,ni->i",
                perturbation_centered[:, :width],
                perturbation_centered[:, offset:],
                optimize=True,
            ) / (valid_pairs - 1)
            denom = np.sqrt(np.clip(variance[:width] * variance[offset:], 1e-24, None))
            corr_band = np.clip(np.nan_to_num(cov_band / denom, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
            diagonal_band_covariance[center + offset, :width] = cov_band.astype(np.float32)
            diagonal_band_correlation[center + offset, :width] = corr_band.astype(np.float32)
            diagonal_band_covariance[center - offset, offset:] = cov_band.astype(np.float32)
            diagonal_band_correlation[center - offset, offset:] = corr_band.astype(np.float32)
    else:
        band_offsets = np.empty((0,), dtype=np.int32)
        diagonal_band_covariance = np.empty((0, 0), dtype=np.float32)
        diagonal_band_correlation = np.empty((0, 0), dtype=np.float32)

    diagonal_block_size = min(max(int(diagonal_block_size), 1), feature_dim)
    if compute_diagonal_blocks:
        diagonal_block_starts_arr = _make_diagonal_block_starts(
            feature_dim,
            diagonal_block_size,
            int(diagonal_block_count),
            diagonal_block_starts,
        )
        diagonal_block_covariance = np.empty(
            (diagonal_block_starts_arr.size, diagonal_block_size, diagonal_block_size),
            dtype=np.float32,
        )
        diagonal_block_correlation = np.empty_like(diagonal_block_covariance)
        for block_index, start in enumerate(tqdm(diagonal_block_starts_arr, desc="Diagonal block statistics", unit="block")):
            block = perturbation_centered[:, int(start) : int(start) + diagonal_block_size]
            cov_block = (block.T @ block) / (valid_pairs - 1)
            cov_block = np.nan_to_num(cov_block, nan=0.0, posinf=0.0, neginf=0.0)
            std = np.sqrt(np.clip(np.diag(cov_block), 1e-12, None))
            corr_block = cov_block / np.outer(std, std)
            corr_block = np.clip(np.nan_to_num(corr_block, nan=0.0, posinf=0.0, neginf=0.0), -1.0, 1.0)
            diagonal_block_covariance[block_index] = cov_block.astype(np.float32)
            diagonal_block_correlation[block_index] = corr_block.astype(np.float32)
    else:
        diagonal_block_starts_arr = np.empty((0,), dtype=np.int32)
        diagonal_block_covariance = np.empty((0, diagonal_block_size, diagonal_block_size), dtype=np.float32)
        diagonal_block_correlation = np.empty_like(diagonal_block_covariance)

    if channel_spatial_diagnostics:
        channel_correlation, spatial_correlation, spatial_correlations, normalized_channel_axis, selected_channel = _compute_channel_spatial_correlations(
            perturbation_centered=perturbation_centered,
            latent_shape=tuple(latent_shape),
            channel_axis=latent_channel_axis,
            spatial_channel_index=spatial_channel_index,
            all_spatial_channel_correlations=all_spatial_channel_correlations,
        )
    else:
        channel_correlation = np.empty((0, 0), dtype=np.float32)
        spatial_correlation = np.empty((0, 0), dtype=np.float32)
        spatial_correlations = np.empty((0, 0, 0), dtype=np.float32)
        normalized_channel_axis = int(latent_channel_axis)
        selected_channel = int(spatial_channel_index)

    return {
        "latent_mean": latent_matrix.mean(axis=0),
        "perturbation_mean": perturbation_mean,
        "variance": variance,
        "B_diag": variance,
        "background_covariance_diag": variance,
        "local_covariance": local_covariance,
        "local_correlation": local_correlation,
        "local_eigenvalues": local_eigenvalues,
        "local_indices": local_indices,
        "latent_shape": np.array(latent_shape, dtype=np.int32),
        "perturbation_centered": perturbation_centered,
        "diagonal_band_covariance": diagonal_band_covariance,
        "diagonal_band_correlation": diagonal_band_correlation,
        "diagonal_band_offsets": band_offsets,
        "diagonal_block_covariance": diagonal_block_covariance,
        "diagonal_block_correlation": diagonal_block_correlation,
        "diagonal_block_starts": diagonal_block_starts_arr,
        "channel_correlation": channel_correlation,
        "spatial_correlation": spatial_correlation,
        "spatial_correlations": spatial_correlations,
        "latent_channel_axis": normalized_channel_axis,
        "spatial_channel_index": selected_channel,
        "sample_count": latent_matrix.shape[0],
        "valid_pairs": valid_pairs,
        "latent_dim": feature_dim,
        "expected_step_hours": expected_step.total_seconds() / 3600.0,
        "expected_pair_step_hours": expected_pair_step.total_seconds() / 3600.0,
    }


def plot_local_covariance(local_covariance: np.ndarray, local_indices: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 6.5))
    cov_limit = np.percentile(np.abs(local_covariance), 99)
    cov_limit = cov_limit if cov_limit > 0 else 1e-6
    cov_img = ax.imshow(local_covariance, cmap="coolwarm", vmin=-cov_limit, vmax=cov_limit, aspect="auto")
    ax.set_title("Local Diagonal Block of Bz")
    ax.set_xlabel(f"Flattened dim ({local_indices[0]}:{local_indices[-1]})")
    ax.set_ylabel(f"Flattened dim ({local_indices[0]}:{local_indices[-1]})")
    fig.colorbar(cov_img, ax=ax, fraction=0.075, pad=0.02, aspect=12)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_correlation_matrix(
    correlation: np.ndarray,
    save_path: Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    correlation = np.asarray(correlation, dtype=np.float32)
    if correlation.size == 0:
        return
    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    image = ax.imshow(
        correlation,
        cmap="bwr",
        vmin=-1.0,
        vmax=1.0,
        origin="upper",
        aspect="equal",
        interpolation="nearest",
        rasterized=True,
    )
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label("Correlation", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    fig.savefig(save_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_paper_correlation_summary(
    channel_correlation: np.ndarray,
    spatial_correlation: np.ndarray,
    save_path: Path,
    *,
    spatial_channel_index: int,
) -> None:
    channel_correlation = np.asarray(channel_correlation, dtype=np.float32)
    spatial_correlation = np.asarray(spatial_correlation, dtype=np.float32)
    if channel_correlation.size == 0 or spatial_correlation.size == 0:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0), constrained_layout=True)
    variable_image = axes[0].imshow(
        channel_correlation,
        cmap="bwr",
        vmin=-1.0,
        vmax=1.0,
        origin="upper",
        aspect="equal",
        interpolation="nearest",
        rasterized=True,
    )
    axes[0].set_title("Variable correlation", fontsize=13, weight="bold")
    axes[0].set_xlabel("Latent channel", fontsize=10)
    axes[0].set_ylabel("Latent channel", fontsize=10)

    axes[1].imshow(
        spatial_correlation,
        cmap="bwr",
        vmin=-1.0,
        vmax=1.0,
        origin="upper",
        aspect="equal",
        interpolation="nearest",
        rasterized=True,
    )
    axes[1].set_title(f"Spatial correlation (latent channel {int(spatial_channel_index)})", fontsize=13, weight="bold")
    axes[1].set_xlabel("Flattened spatial index", fontsize=10)
    axes[1].set_ylabel("Flattened spatial index", fontsize=10)

    colorbar = fig.colorbar(variable_image, ax=axes, fraction=0.032, pad=0.018)
    colorbar.set_label("Correlation", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    fig.savefig(save_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def _latent_shape_tuple(latent_shape: np.ndarray | tuple[int, ...]) -> tuple[int, ...]:
    if isinstance(latent_shape, np.ndarray):
        return tuple(int(v) for v in latent_shape.tolist())
    return tuple(int(v) for v in latent_shape)


def _reshape_flat_field_for_plot(values: np.ndarray, latent_shape: np.ndarray | tuple[int, ...]) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    shape = _latent_shape_tuple(latent_shape)
    if len(shape) > 1 and int(np.prod(shape)) == values.size:
        field = values.reshape(shape)
        while field.ndim > 2:
            field = field.mean(axis=-1)
        return field

    rows = int(np.floor(np.sqrt(values.size)))
    rows = max(rows, 1)
    cols = int(np.ceil(values.size / rows))
    padded = np.full(rows * cols, np.nan, dtype=np.float64)
    padded[: values.size] = values
    return padded.reshape(rows, cols)


def plot_bz_full_field(variance: np.ndarray, latent_shape: np.ndarray | tuple[int, ...], save_path: Path) -> None:
    field = _reshape_flat_field_for_plot(variance, latent_shape)
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    image = ax.imshow(field, cmap="magma", aspect="auto")
    ax.set_title("Full-Field Bz Diagonal")
    ax.set_xlabel("Latent x / flattened block")
    ax.set_ylabel("Latent y / flattened block")
    fig.colorbar(image, ax=ax, fraction=0.075, pad=0.02, aspect=12, label="Bz diagonal variance")
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_local_correlation(local_correlation: np.ndarray, local_indices: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 6.5))
    corr_img = ax.imshow(local_correlation, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="auto")
    ax.set_title("Flattened Latent Local Correlation")
    ax.set_xlabel(f"Flattened dim ({local_indices[0]}:{local_indices[-1]})")
    ax.set_ylabel(f"Flattened dim ({local_indices[0]}:{local_indices[-1]})")
    fig.colorbar(corr_img, ax=ax, fraction=0.075, pad=0.02, aspect=12)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_local_eigenvalues(local_eigenvalues: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    top_k = min(100, local_eigenvalues.shape[0])
    eigvals_to_plot = np.clip(local_eigenvalues[:top_k], a_min=1e-12, a_max=None)
    ax.semilogy(np.arange(1, top_k + 1), eigvals_to_plot, linewidth=1.6)
    ax.set_title("Leading Eigenvalues of Local Block")
    ax.set_xlabel("Mode index")
    ax.set_ylabel("Eigenvalue")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_variance_profile(variance: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(variance, linewidth=1.0)
    ax.set_title("Background Variance by Flattened Latent Dimension")
    ax.set_xlabel("Flattened latent dimension")
    ax.set_ylabel("Variance")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_variance_histogram(variance: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(variance, bins=40, color="#2f6c8f", edgecolor="white")
    ax.set_title("Histogram of Flattened Latent Background Variance")
    ax.set_xlabel("Variance")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_variance_spatial_map(variance: np.ndarray, latent_shape: np.ndarray, save_path: Path) -> None:
    latent_shape_tuple = tuple(int(v) for v in latent_shape.tolist())
    if len(latent_shape_tuple) <= 1:
        return

    variance_map = variance.reshape(latent_shape_tuple)
    spatial_map = variance_map
    while spatial_map.ndim > 2:
        spatial_map = spatial_map.mean(axis=-1)

    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    spatial_img = ax.imshow(spatial_map, cmap="magma", aspect="auto")
    ax.set_title("Variance Map Reshaped from Flattened Latent")
    ax.set_xlabel("Latent x")
    ax.set_ylabel("Latent y")
    fig.colorbar(spatial_img, ax=ax, fraction=0.075, pad=0.02, aspect=12)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_variance_channel_profile(variance: np.ndarray, latent_shape: np.ndarray, save_path: Path) -> None:
    latent_shape_tuple = tuple(int(v) for v in latent_shape.tolist())
    if len(latent_shape_tuple) <= 1:
        return

    variance_map = variance.reshape(latent_shape_tuple)
    if variance_map.ndim >= 3:
        reduce_axes = tuple(range(variance_map.ndim - 1))
        channel_profile = variance_map.mean(axis=reduce_axes)
    else:
        channel_profile = variance_map.mean(axis=0)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.plot(channel_profile, linewidth=1.0)
    ax.set_title("Mean Variance by Latent Channel")
    ax.set_xlabel("Latent channel")
    ax.set_ylabel("Variance")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
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
    centered: np.ndarray,
    corr_bins: int,
    diagnostic_pair_count: int,
    pair_batch_size: int,
    random_seed: int,
) -> dict[str, np.ndarray | float | int | str]:
    if centered.ndim != 2:
        raise ValueError("centered perturbations must have shape (N, D).")
    sample_count, feature_dim = centered.shape
    if sample_count < 2 or feature_dim < 2:
        raise RuntimeError("At least two perturbations and two latent dimensions are required.")

    corr_bins = max(3, int(corr_bins))
    pair_batch_size = max(1, int(pair_batch_size))
    signed_edges = np.linspace(-1.0, 1.0, corr_bins + 1, dtype=np.float64)
    abs_edges = np.linspace(0.0, 1.0, corr_bins + 1, dtype=np.float64)
    signed_hist = np.zeros(corr_bins, dtype=np.int64)
    abs_hist = np.zeros(corr_bins, dtype=np.int64)

    sum_sq = np.square(centered, dtype=np.float64).sum(axis=0, dtype=np.float64)
    sum_sq = np.clip(sum_sq, a_min=1e-24, a_max=None)

    def add_values(values: np.ndarray) -> None:
        values = np.clip(np.asarray(values, dtype=np.float64).reshape(-1), -1.0, 1.0)
        signed_hist[:] += np.histogram(values, bins=signed_edges)[0]
        abs_hist[:] += np.histogram(np.abs(values), bins=abs_edges)[0]

    pair_i, pair_j, exact_pairs = _make_offdiag_pairs(feature_dim, int(diagnostic_pair_count), int(random_seed))
    pair_count = int(pair_i.size)
    processed_pairs = 0
    with tqdm(total=pair_count, desc="Offdiag correlation pairs", unit="pair") as progress:
        while processed_pairs < pair_count:
            stop = min(processed_pairs + pair_batch_size, pair_count)
            left = pair_i[processed_pairs:stop]
            right = pair_j[processed_pairs:stop]
            numerator = (centered[:, left] * centered[:, right]).sum(axis=0, dtype=np.float64)
            denom = np.sqrt(sum_sq[left] * sum_sq[right])
            add_values(numerator / denom)
            progress.update(stop - processed_pairs)
            processed_pairs = stop

    mode = "global_all_nmc_pairs" if exact_pairs else "global_random_nmc_pairs"

    signed_p95 = percentile_from_histogram(signed_hist, signed_edges, 95.0)
    abs_p95 = percentile_from_histogram(abs_hist, abs_edges, 95.0)
    return {
        "offdiag_corr_hist": signed_hist,
        "offdiag_corr_bin_edges": signed_edges,
        "offdiag_abs_corr_hist": abs_hist,
        "offdiag_abs_corr_bin_edges": abs_edges,
        "offdiag_corr_p95": signed_p95,
        "offdiag_abs_corr_p95": abs_p95,
        "offdiag_corr_pair_count": int(pair_count),
        "offdiag_corr_distribution_mode": mode,
    }


def plot_offdiag_correlation_distribution(
    corr_distribution: dict[str, np.ndarray | float | int | str],
    save_path: Path,
) -> None:
    hist = np.asarray(corr_distribution["offdiag_corr_hist"], dtype=np.float64)
    bin_edges = np.asarray(corr_distribution["offdiag_corr_bin_edges"], dtype=np.float64)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    density = hist / max(float(hist.sum()), 1.0) / widths
    abs_p95 = float(corr_distribution["offdiag_abs_corr_p95"])
    mode = str(corr_distribution["offdiag_corr_distribution_mode"])

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(centers, density, color="#1b6ca8", linewidth=2.0)
    ax.axvline(0.0, color="#333333", linewidth=1.2, alpha=0.7)
    ax.axvline(abs_p95, color="#f18f01", linewidth=2.0, linestyle="--", label=f"95% |corr| = {abs_p95:.3f}")
    ax.axvline(-abs_p95, color="#f18f01", linewidth=2.0, linestyle="--")
    ax.set_xlim(-1.0, 1.0)
    ax.set_title(f"Off-Diagonal Correlation Distribution ({mode})")
    ax.set_xlabel("corr(i, j), i != j")
    ax.set_ylabel("Density")
    ax.grid(True, alpha=0.25)
    ax.legend()
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _robust_symmetric_limit(values: np.ndarray, percentile: float = 99.7) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 1.0
    limit = float(np.percentile(np.abs(finite), percentile))
    if not np.isfinite(limit) or limit <= 0.0:
        limit = float(np.max(np.abs(finite))) if finite.size else 1.0
    return max(limit, 1e-12)


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
    save_path: Path,
    *,
    title: str,
    label: str,
    tile_width: int = 1024,
    cmap: str = "bwr",
) -> None:
    band = np.asarray(band, dtype=np.float32)
    offsets = np.asarray(offsets, dtype=np.int32)
    mosaic = _make_band_mosaic(band, tile_width=tile_width)
    finite = mosaic[np.isfinite(mosaic)]
    limit = _robust_symmetric_limit(finite, percentile=99.7) if finite.size else 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
    fig_height = float(np.clip(mosaic.shape[0] / 85.0, 5.0, 24.0))
    fig_width = float(np.clip(int(tile_width) / 260.0, 8.0, 18.0))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    image = ax.imshow(
        mosaic,
        origin="upper",
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        rasterized=True,
    )
    rows = band.shape[0]
    n_tiles = int(np.ceil(band.shape[1] / max(int(tile_width), 1)))
    centers = [tile * (rows + 2) + rows / 2 - 0.5 for tile in range(n_tiles)]
    labels = [f"{tile * int(tile_width):,}" for tile in range(n_tiles)]
    if len(centers) > 18:
        step = int(np.ceil(len(centers) / 18))
        centers = centers[::step]
        labels = labels[::step]
    ax.set_yticks(centers)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"Tile-local latent index, tile width={int(tile_width):,}", fontsize=10)
    ax.set_ylabel("Tile start flattened latent index", fontsize=10)
    ax.set_title(f"{title}\noffsets {int(offsets[0])}..{int(offsets[-1])}", fontsize=12, pad=9)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.018)
    colorbar.set_label(label, fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    colorbar.formatter.set_powerlimits((-2, 2))
    colorbar.update_ticks()
    fig.savefig(save_path, dpi=PAPER_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_diagonal_block_mosaic(
    blocks: np.ndarray,
    starts: np.ndarray,
    save_path: Path,
    *,
    title: str,
    label: str,
    tile_cols: int = 4,
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
    fig.savefig(save_path, dpi=PAPER_DPI, bbox_inches="tight")
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
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
        colorbar.set_label(label, fontsize=10)
        colorbar.ax.tick_params(labelsize=8)
        colorbar.formatter.set_powerlimits((-2, 2))
        colorbar.update_ticks()
        fig.savefig(output_dir / f"{prefix}_{block_index:03d}_start_{int(start):06d}.png", dpi=PAPER_DPI, bbox_inches="tight")
        plt.close(fig)


def compute_diag_approx_diagnostics(local_covariance: np.ndarray, local_correlation: np.ndarray) -> dict[str, np.ndarray | float]:
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

    row_offdiag_abs_sum = np.abs(offdiag_entries).sum(axis=1)
    row_coupling_ratio = row_offdiag_abs_sum / (np.abs(diag_entries) + 1e-12)

    return {
        "offdiag_abs_corr": offdiag_abs_corr,
        "diag_energy": diag_energy,
        "offdiag_energy": offdiag_energy,
        "diag_energy_fraction": float(diag_energy / total_energy) if total_energy > 0 else 0.0,
        "offdiag_energy_fraction": float(offdiag_energy / total_energy) if total_energy > 0 else 0.0,
        "row_coupling_ratio": row_coupling_ratio,
        "offdiag_abs_corr_median": float(np.median(offdiag_abs_corr)) if offdiag_abs_corr.size else 0.0,
        "offdiag_abs_corr_p95": float(np.percentile(offdiag_abs_corr, 95)) if offdiag_abs_corr.size else 0.0,
        "row_coupling_ratio_median": float(np.median(row_coupling_ratio)) if row_coupling_ratio.size else 0.0,
        "row_coupling_ratio_p95": float(np.percentile(row_coupling_ratio, 95)) if row_coupling_ratio.size else 0.0,
    }


def plot_offdiag_corr_histogram(offdiag_abs_corr: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(offdiag_abs_corr, bins=50, color="#3f7f93", edgecolor="white")
    if offdiag_abs_corr.size:
        median = float(np.median(offdiag_abs_corr))
        p95 = float(np.percentile(offdiag_abs_corr, 95))
        ax.axvline(median, color="#a23b72", linewidth=2, label=f"median={median:.3g}")
        ax.axvline(p95, color="#f18f01", linewidth=2, linestyle="--", label=f"p95={p95:.3g}")
    ax.set_title("Absolute Off-Diagonal Correlation Histogram")
    ax.set_xlabel("|corr(i, j)|, i != j")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    if offdiag_abs_corr.size:
        ax.legend()
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_offdiag_corr_cdf(offdiag_abs_corr: np.ndarray, save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if offdiag_abs_corr.size:
        sorted_vals = np.sort(offdiag_abs_corr)
        cdf = np.linspace(0.0, 1.0, sorted_vals.size, endpoint=True)
        ax.plot(sorted_vals, cdf, linewidth=2.0, color="#1b6ca8")
    ax.set_title("CDF of Absolute Off-Diagonal Correlation")
    ax.set_xlabel("|corr(i, j)|, i != j")
    ax.set_ylabel("CDF")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_diag_energy_fraction(diag_stats: dict[str, np.ndarray | float], save_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    diag_energy_fraction = float(diag_stats["diag_energy_fraction"])
    offdiag_energy_fraction = float(diag_stats["offdiag_energy_fraction"])
    ax.bar(["diag", "offdiag"], [diag_energy_fraction, offdiag_energy_fraction], color=["#2a9d8f", "#e76f51"])
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Frobenius Energy Fraction in Local Block")
    ax.set_ylabel("Energy fraction")
    for idx, value in enumerate([diag_energy_fraction, offdiag_energy_fraction]):
        ax.text(idx, value + 0.02, f"{value:.3f}", ha="center", va="bottom")
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_row_coupling_ratio_histogram(diag_stats: dict[str, np.ndarray | float], save_path: Path) -> None:
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
    plt.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.lag < 1:
        raise ValueError("--lag must be >= 1.")
    if args.plot_dims < 1:
        raise ValueError("--plot-dims must be >= 1.")
    if args.corr_bins < 3:
        raise ValueError("--corr-bins must be >= 3.")
    if args.corr_sample_pairs is not None:
        args.diagnostic_pair_count = int(args.corr_sample_pairs)
    if args.corr_pair_batch_size is not None:
        args.diagnostic_pair_batch_size = int(args.corr_pair_batch_size)
    if args.corr_random_seed is not None:
        args.diagnostic_seed = int(args.corr_random_seed)
    if args.diagnostic_pair_count < 1:
        raise ValueError("--diagnostic-pair-count must be >= 1.")
    if args.diagnostic_pair_batch_size < 1:
        raise ValueError("--diagnostic-pair-batch-size must be >= 1.")
    if args.diagonal_block_size < 1:
        raise ValueError("--diagonal-block-size must be >= 1.")
    if args.diagonal_block_count < 0:
        raise ValueError("--diagonal-block-count must be >= 0.")
    if args.plot_channel_spatial_correlation and not args.channel_spatial_diagnostics:
        raise ValueError("--plot-channel-spatial-correlation requires --channel-spatial-diagnostics.")
    if args.plot_all_spatial_channel_correlations and not args.channel_spatial_diagnostics:
        raise ValueError("--plot-all-spatial-channel-correlations requires --channel-spatial-diagnostics.")

    device = torch.device(args.device or ("cuda" if torch.cuda.is_available() else "cpu"))

    dataset = WRFSlidingDataset(
        nc_dir=args.nc_dir,
        input_len=DATA_CONFIG["input_len"],
        label_len=DATA_CONFIG["label_len"],
        input_vars=DATA_CONFIG["input_vars"],
        label_vars=DATA_CONFIG["label_vars"],
    )
    stats_subset, stats_indices, sample_times, selected_quarters = build_stats_subset(
        dataset=dataset,
        base_year=args.base_year,
        quarters=args.quarters,
        max_samples=args.max_samples,
    )

    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else (
        Path(__file__).resolve().parent / f"nmc_outputs_q{'_'.join(map(str, selected_quarters))}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    dataloader = DataLoader(
        stats_subset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=(device.type == "cuda"),
    )

    model_path = resolve_model_path(args.model_path)
    model = build_model().to(device)
    load_checkpoint(model, model_path, device)
    model.eval()

    print(f"[INFO] model checkpoint: {model_path}")
    print(f"[INFO] statistics quarters: {selected_quarters}")
    print(f"[INFO] statistics samples: {len(stats_subset)}")
    print(f"[INFO] first selected sample index: {stats_indices[0]}")
    print(f"[INFO] last selected sample index: {stats_indices[-1]}")

    stats = collect_nmc_statistics(
        model=model,
        dataloader=dataloader,
        sample_times=sample_times,
        device=device,
        lag=args.lag,
        plot_dims=args.plot_dims,
        plot_offset=args.plot_offset,
        diagonal_band_radius=args.diagonal_band_radius,
        diagonal_block_size=args.diagonal_block_size,
        diagonal_block_count=args.diagonal_block_count,
        diagonal_block_starts=args.diagonal_block_starts,
        compute_diagonal_band=args.plot_diagonal_band,
        compute_diagonal_blocks=args.plot_diagonal_blocks or args.save_diagonal_block_tiles,
        channel_spatial_diagnostics=args.channel_spatial_diagnostics,
        latent_channel_axis=args.latent_channel_axis,
        spatial_channel_index=args.spatial_channel_index,
        all_spatial_channel_correlations=args.plot_all_spatial_channel_correlations,
    )

    bz_full_field_path = output_dir / "latent_nmc_Bz_full_field.png"
    bz_local_block_path = output_dir / "latent_nmc_Bz_local_diagonal_block.png"
    offdiag_corr_distribution_path = output_dir / "latent_nmc_offdiag_correlation_distribution.png"
    bz_diagonal_band_path = output_dir / "latent_nmc_Bz_diagonal_band_mosaic.png"
    corr_diagonal_band_path = output_dir / "latent_nmc_corr_diagonal_band_mosaic.png"
    bz_diagonal_block_path = output_dir / "latent_nmc_Bz_diagonal_block_mosaic.png"
    corr_diagonal_block_path = output_dir / "latent_nmc_corr_diagonal_block_mosaic.png"
    channel_correlation_path = output_dir / "latent_nmc_channel_correlation.png"
    spatial_correlation_path = output_dir / "latent_nmc_spatial_correlation.png"
    paper_correlation_summary_path = output_dir / "latent_nmc_paper_correlation_summary.png"
    all_spatial_correlation_dir = output_dir / "latent_nmc_spatial_correlations_all"

    plot_bz_full_field(stats["B_diag"], stats["latent_shape"], bz_full_field_path)
    if args.plot_local_block:
        plot_local_covariance(stats["local_covariance"], stats["local_indices"], bz_local_block_path)
    diag_stats = compute_diag_approx_diagnostics(stats["local_covariance"], stats["local_correlation"])
    corr_distribution = compute_offdiag_correlation_distribution(
        centered=stats["perturbation_centered"],
        corr_bins=args.corr_bins,
        diagnostic_pair_count=args.diagnostic_pair_count,
        pair_batch_size=args.diagnostic_pair_batch_size,
        random_seed=args.diagnostic_seed,
    )
    plot_offdiag_correlation_distribution(corr_distribution, offdiag_corr_distribution_path)
    if args.plot_channel_spatial_correlation:
        plot_correlation_matrix(
            stats["channel_correlation"],
            channel_correlation_path,
            title="Latent Channel Correlation",
            xlabel="Latent channel",
            ylabel="Latent channel",
        )
        plot_correlation_matrix(
            stats["spatial_correlation"],
            spatial_correlation_path,
            title=f"Spatial Correlation (latent channel {int(stats['spatial_channel_index'])})",
            xlabel="Flattened spatial index",
            ylabel="Flattened spatial index",
        )
        plot_paper_correlation_summary(
            stats["channel_correlation"],
            stats["spatial_correlation"],
            paper_correlation_summary_path,
            spatial_channel_index=int(stats["spatial_channel_index"]),
        )
    if args.plot_all_spatial_channel_correlations:
        spatial_correlations = np.asarray(stats["spatial_correlations"], dtype=np.float32)
        all_spatial_correlation_dir.mkdir(parents=True, exist_ok=True)
        for channel_index in tqdm(range(spatial_correlations.shape[0]), desc="Plotting spatial channel correlations", unit="channel"):
            plot_correlation_matrix(
                spatial_correlations[channel_index],
                all_spatial_correlation_dir / f"latent_nmc_spatial_correlation_channel_{channel_index:03d}.png",
                title=f"Spatial Correlation (latent channel {channel_index})",
                xlabel="Flattened spatial index",
                ylabel="Flattened spatial index",
            )
    if args.plot_diagonal_band:
        plot_diagonal_band_mosaic(
            stats["diagonal_band_covariance"],
            stats["diagonal_band_offsets"],
            bz_diagonal_band_path,
            title="Full diagonal band mosaic of NMC Bz",
            label="Covariance",
            tile_width=args.diagonal_band_tile_width,
            cmap="bwr",
        )
        plot_diagonal_band_mosaic(
            stats["diagonal_band_correlation"],
            stats["diagonal_band_offsets"],
            corr_diagonal_band_path,
            title="Full diagonal band mosaic of NMC correlation",
            label="Correlation",
            tile_width=args.diagonal_band_tile_width,
            cmap="bwr",
        )
    if args.plot_diagonal_blocks:
        plot_diagonal_block_mosaic(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            bz_diagonal_block_path,
            title=f"Diagonal {args.diagonal_block_size:,}x{args.diagonal_block_size:,} blocks of NMC Bz",
            label="Covariance",
            tile_cols=args.diagonal_block_tile_cols,
            cmap="bwr",
        )
        plot_diagonal_block_mosaic(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            corr_diagonal_block_path,
            title=f"Diagonal {args.diagonal_block_size:,}x{args.diagonal_block_size:,} blocks of NMC correlation",
            label="Correlation",
            tile_cols=args.diagonal_block_tile_cols,
            cmap="bwr",
        )
    if args.save_diagonal_block_tiles:
        block_dir = output_dir / "diagonal_block_tiles"
        plot_diagonal_block_tiles(
            stats["diagonal_block_covariance"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="Bz_covariance",
            title="NMC Bz diagonal block",
            label="Covariance",
            cmap="bwr",
        )
        plot_diagonal_block_tiles(
            stats["diagonal_block_correlation"],
            stats["diagonal_block_starts"],
            block_dir,
            prefix="correlation",
            title="NMC correlation diagonal block",
            label="Correlation",
            cmap="bwr",
        )

    npz_path = output_dir / "latent_nmc_background_covariance.npz"
    np.savez_compressed(
        npz_path,
        # Complete diagonal NMC data used by WRF-VAE-LDA. The local block below is
        # only a diagnostic view around the diagonal.
        variance=stats["variance"],
        B_diag=stats["B_diag"],
        background_covariance_diag=stats["background_covariance_diag"],
        local_covariance=stats["local_covariance"],
        local_correlation=stats["local_correlation"],
        local_eigenvalues=stats["local_eigenvalues"],
        local_indices=stats["local_indices"],
        latent_mean=stats["latent_mean"],
        perturbation_mean=stats["perturbation_mean"],
        latent_shape=stats["latent_shape"],
        diagonal_band_radius=np.array([args.diagonal_band_radius], dtype=np.int32),
        diagonal_block_size=np.array([args.diagonal_block_size], dtype=np.int32),
        diagonal_block_starts=stats["diagonal_block_starts"],
        channel_correlation=stats["channel_correlation"],
        spatial_correlation=stats["spatial_correlation"],
        spatial_correlations=stats["spatial_correlations"],
        latent_channel_axis=np.array([stats["latent_channel_axis"]], dtype=np.int32),
        spatial_channel_index=np.array([stats["spatial_channel_index"]], dtype=np.int32),
        local_offdiag_abs_corr=diag_stats["offdiag_abs_corr"],
        local_row_coupling_ratio=diag_stats["row_coupling_ratio"],
        offdiag_corr_hist=corr_distribution["offdiag_corr_hist"],
        offdiag_corr_bin_edges=corr_distribution["offdiag_corr_bin_edges"],
        offdiag_abs_corr_hist=corr_distribution["offdiag_abs_corr_hist"],
        offdiag_abs_corr_bin_edges=corr_distribution["offdiag_abs_corr_bin_edges"],
        offdiag_corr_p95=np.array([corr_distribution["offdiag_corr_p95"]], dtype=np.float64),
        offdiag_abs_corr_p95=np.array([corr_distribution["offdiag_abs_corr_p95"]], dtype=np.float64),
        offdiag_corr_pair_count=np.array([corr_distribution["offdiag_corr_pair_count"]], dtype=np.int64),
        offdiag_corr_distribution_mode=np.array([corr_distribution["offdiag_corr_distribution_mode"]], dtype=object),
        diagnostic_pair_count=np.array([corr_distribution["offdiag_corr_pair_count"]], dtype=np.int64),
        diagnostic_seed=np.array([args.diagnostic_seed], dtype=np.int64),
        diag_energy_fraction=np.array([diag_stats["diag_energy_fraction"]], dtype=np.float64),
        offdiag_energy_fraction=np.array([diag_stats["offdiag_energy_fraction"]], dtype=np.float64),
        offdiag_abs_corr_median=np.array([diag_stats["offdiag_abs_corr_median"]], dtype=np.float64),
        local_offdiag_abs_corr_p95=np.array([diag_stats["offdiag_abs_corr_p95"]], dtype=np.float64),
        row_coupling_ratio_median=np.array([diag_stats["row_coupling_ratio_median"]], dtype=np.float64),
        row_coupling_ratio_p95=np.array([diag_stats["row_coupling_ratio_p95"]], dtype=np.float64),
        sample_count=np.array([stats["sample_count"]], dtype=np.int32),
        valid_pairs=np.array([stats["valid_pairs"]], dtype=np.int32),
        latent_dim=np.array([stats["latent_dim"]], dtype=np.int32),
        expected_step_hours=np.array([stats["expected_step_hours"]], dtype=np.float32),
        expected_pair_step_hours=np.array([stats["expected_pair_step_hours"]], dtype=np.float32),
        selected_quarters=np.array(selected_quarters, dtype=np.int32),
    )

    print(f"[INFO] latent dim: {stats['latent_dim']}")
    print(f"[INFO] valid NMC pairs: {stats['valid_pairs']}")
    print(f"[INFO] local covariance block: [{stats['local_indices'][0]}:{stats['local_indices'][-1]}]")
    print(f"[INFO] diag energy fraction (local block): {diag_stats['diag_energy_fraction']:.4f}")
    print(f"[INFO] offdiag |corr| p95 ({corr_distribution['offdiag_corr_distribution_mode']} distribution): {corr_distribution['offdiag_abs_corr_p95']:.4f}")
    print(f"[INFO] inferred sample step: {stats['expected_step_hours']:.3f} hours")
    print(f"[INFO] inferred pair step: {stats['expected_pair_step_hours']:.3f} hours")
    print(f"[INFO] covariance file: {npz_path}")
    print(f"[INFO] Bz full-field figure: {bz_full_field_path}")
    if args.plot_local_block:
        print(f"[INFO] Bz local diagonal block figure: {bz_local_block_path}")
    print(f"[INFO] offdiag correlation distribution figure: {offdiag_corr_distribution_path}")
    if args.plot_channel_spatial_correlation:
        print(f"[INFO] channel correlation figure: {channel_correlation_path}")
        print(f"[INFO] spatial correlation figure: {spatial_correlation_path}")
        print(f"[INFO] paper correlation summary figure: {paper_correlation_summary_path}")
    if args.plot_all_spatial_channel_correlations:
        print(f"[INFO] all spatial channel correlation figures: {all_spatial_correlation_dir}")
    if args.plot_diagonal_band:
        print(f"[INFO] Bz diagonal band mosaic figure: {bz_diagonal_band_path}")
    if args.plot_diagonal_blocks:
        print(f"[INFO] Bz diagonal block mosaic figure: {bz_diagonal_block_path}")
    if args.save_diagonal_block_tiles:
        print(f"[INFO] diagonal block tile figures: {output_dir / 'diagonal_block_tiles'}")


if __name__ == "__main__":
    main()
