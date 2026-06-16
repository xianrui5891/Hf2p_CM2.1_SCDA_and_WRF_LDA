from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def _field_figsize(shape: tuple[int, int], ncols: int = 1, base_height: float = 4.0) -> tuple[float, float]:
    height, width = int(shape[0]), int(shape[1])
    aspect = width / max(height, 1)
    panel_width = float(np.clip(base_height * aspect, 3.5, 7.2))
    colorbar_padding = 0.7 * ncols
    return panel_width * ncols + colorbar_padding, base_height


def _as_masked_field(data: np.ndarray) -> np.ma.MaskedArray:
    if np.ma.isMaskedArray(data):
        return np.ma.masked_invalid(data)
    return np.ma.masked_invalid(np.asarray(data))


def _field_cmap(name: str):
    cmap = plt.get_cmap(name)
    try:
        cmap = cmap.copy()
    except AttributeError:
        pass
    cmap.set_bad(color="#f2f2f2", alpha=0.0)
    return cmap


def _runtime_cmap(var_name: str, data_type: str, default: str) -> str:
    name = var_name.lower()
    if data_type == "increment":
        return default
    if "residual" in name or "improvement" in name:
        return default
    if name in {"sst", "sea_surface_temperature"}:
        return "turbo"
    return default


def _runtime_limits(var_name: str, data_type: str, values: np.ndarray) -> tuple[float, float]:
    name = var_name.lower()
    if data_type == "increment" or "residual" in name or "improvement" in name:
        bound = max(float(np.nanpercentile(np.abs(values), 98.5)), 1e-12)
        return -bound, bound
    if name in {"sst", "sea_surface_temperature"}:
        lo = float(np.nanpercentile(values, 0.5))
        hi = float(np.nanpercentile(values, 99.5))
    else:
        lo = float(np.nanpercentile(values, 1.0))
        hi = float(np.nanpercentile(values, 99.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        lo = float(np.nanmin(values))
        hi = float(np.nanmax(values))
    return lo, hi


def _finite_abs_bound(*arrays: np.ndarray, fallback: float = 1e-12) -> float:
    values = []
    for arr in arrays:
        masked = _as_masked_field(arr)
        compressed = np.asarray(masked.compressed(), dtype=np.float64)
        if compressed.size:
            values.append(np.abs(compressed))
    if not values:
        return fallback
    bound = float(np.max(np.concatenate(values)))
    return max(bound, fallback)


class FieldPlotter:
    """Field plotting compatible with temp/cm2_cda_main.py."""

    def __init__(self, base_dir: str | Path, cmap: str = "bwr", use_runtime_cmap: bool = True) -> None:
        self.base_dir = Path(base_dir)
        self.cmap = cmap
        self.use_runtime_cmap = bool(use_runtime_cmap)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def plot_and_save(
        self,
        time_str: str,
        var_name: str,
        data2d: np.ndarray,
        data_type: str = "original",
        vmin=None,
        vmax=None,
    ) -> str:
        folder = self.base_dir / data_type / var_name
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / f"{time_str}_{data_type}_{var_name}.png"
        arr = _as_masked_field(data2d)
        if arr.ndim != 2:
            raise ValueError("data2d must be 2D")
        display = arr.T
        fig, ax = plt.subplots(figsize=_field_figsize(display.shape, base_height=5.0))
        if vmin is None and vmax is None and display.count() > 0:
            values = np.asarray(display.compressed(), dtype=np.float64)
            vmin, vmax = _runtime_limits(var_name, data_type, values)
        cmap_name = _runtime_cmap(var_name, data_type, self.cmap) if self.use_runtime_cmap else self.cmap
        image = ax.imshow(display, origin="lower", aspect="equal", cmap=_field_cmap(cmap_name), vmin=vmin, vmax=vmax)
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_title(f"{time_str} {data_type} {var_name}")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.12)
        fig.colorbar(image, cax=cax, label=var_name)
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return str(path)


class LossPlotter:
    def __init__(self, output_dir: str | Path) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot(self, history: Sequence[Mapping[str, float]], filename: str = "losses.png") -> None:
        if not history:
            return
        self.output_dir.mkdir(parents=True, exist_ok=True)

        def series(key: str) -> list[float]:
            return [float(item[key]) for item in history if key in item and np.isfinite(float(item[key]))]

        metric_names = {
            key.removeprefix("train_").removeprefix("val_")
            for row in history
            for key in row
            if (key.startswith("train_") or key.startswith("val_"))
            and (
                key.endswith("_weighted_loss")
                or key.endswith("_total")
                or key.endswith("total_loss")
            )
            and not key.removeprefix("train_").removeprefix("val_").startswith("corrupt_")
        }
        preferred = [
            "total_loss",
            "atm_total",
            "ocn_total",
            "atm_l1_weighted_loss",
            "atm_l2_weighted_loss",
            "atm_subband_weighted_loss",
            "ocn_l1_weighted_loss",
            "ocn_l2_weighted_loss",
            "ocn_subband_weighted_loss",
            "atm_spectral_weighted_loss",
            "ocn_spectral_weighted_loss",
            "atm_gradient_weighted_loss",
            "ocn_gradient_weighted_loss",
            "atm_energy_weighted_loss",
            "ocn_energy_weighted_loss",
        ]
        ordered = [key for key in preferred if key in metric_names] + sorted(
            {
                key
                for key in metric_names
                if key not in preferred
            }
        )
        filtered = []
        for key in ordered:
            values = series(f"train_{key}") + series(f"val_{key}")
            if key in {"total_loss", "atm_subband_weighted_loss", "ocn_subband_weighted_loss"} or (
                values and max(abs(v) for v in values) >= 1e-4
            ):
                filtered.append(key)
        keys = filtered[:15]
        if not keys:
            return
        ncols = min(2, len(keys))
        nrows = int(np.ceil(len(keys) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 3.8 * nrows), squeeze=False)
        flat = axes.ravel()
        for ax, key in zip(flat, keys):
            train_values = series(f"train_{key}")
            val_values = series(f"val_{key}")
            if train_values:
                ax.plot(train_values, label="train", color="#1f77b4", linewidth=1.5)
            if val_values:
                ax.plot(val_values, label="val", color="#d62728", linewidth=1.5)
            ax.set_title(key)
            ax.grid(True, alpha=0.3)
            if key == "total_loss":
                ax.set_yscale("log")
            if train_values and val_values:
                ax.legend()
        for ax in flat[len(keys) :]:
            ax.axis("off")
        fig.tight_layout()
        output_path = self.output_dir / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150)
        plt.close(fig)


def plot_reconstruction_triplet(
    original: np.ndarray,
    reconstructed: np.ndarray,
    output_path: str | Path,
    title: str,
    cmap: str = "bwr",
) -> None:
    original = _as_masked_field(original)
    reconstructed = _as_masked_field(reconstructed)
    shared_mask = np.ma.getmaskarray(original) | np.ma.getmaskarray(reconstructed)
    original = np.ma.array(original, mask=shared_mask)
    reconstructed = np.ma.array(reconstructed, mask=shared_mask)
    diff = np.ma.masked_invalid(reconstructed - original)
    vmax = _finite_abs_bound(original, reconstructed)
    err = _finite_abs_bound(diff)
    fig, axes = plt.subplots(1, 3, figsize=_field_figsize(original.shape, ncols=3, base_height=4.0))
    panels = [
        (original, "Original", vmax),
        (reconstructed, "Reconstructed", vmax),
        (diff, "Error", err),
    ]
    image_cmap = _field_cmap(cmap)
    for ax, (arr, label, bound) in zip(axes, panels):
        image = ax.imshow(arr, origin="lower", aspect="equal", cmap=image_cmap, vmin=-bound, vmax=bound)
        ax.set_title(label)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_spectrum_pair(
    wavenumber: np.ndarray,
    original_spectrum: np.ndarray,
    reconstructed_spectrum: np.ndarray,
    output_path: str | Path,
    title: str,
    colors: tuple[str, str] = ("#1f77b4", "#d62728"),
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(wavenumber, original_spectrum, label="Original", linewidth=1.8, color=colors[0])
    ax.plot(wavenumber, reconstructed_spectrum, label="Reconstructed", linewidth=1.8, color=colors[1])
    ax.set_yscale("log")
    ax.set_xlabel("Normalized wavenumber")
    ax.set_ylabel("Radial power")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
