from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset

PROJECT_DIR = Path(__file__).resolve().parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from dataset import WRFSlidingDataset
from utils.normalization_utils import norm_dict, var_name as DEFAULT_VARIABLES
from vaelda.networks.transformer import LGUnet_all

DEFAULT_MODEL_PATH = PROJECT_DIR / "model_forgoal.pth"
DEFAULT_DATA_DIR = Path("/data/zsq_data/d01_d02_d03/d03_wrf")


@dataclass
class VariableMetric:
    domain: str
    variable: str
    sample: int
    mae: float
    median_ae: float
    p95_abs_error: float
    mse: float
    rmse: float
    centered_rmse: float
    nmae_range: float
    nrmse_range: float
    nrmse_std: float
    bias: float
    smape: float
    mape: float
    max_abs_error: float
    relative_l2: float
    correlation: float
    r2: float
    explained_variance: float
    nse: float
    kge: float
    willmott_d: float
    ssim_global: float
    psnr: float
    spectral_l1: float


def radial_power_spectrum_np(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    arr = np.nan_to_num(np.asarray(field, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    height, width = arr.shape
    power = np.abs(np.fft.rfft2(arr, norm="ortho")) ** 2
    ky = np.fft.fftfreq(height)[:, None]
    kx = np.fft.rfftfreq(width)[None, :]
    radius = np.sqrt(ky**2 + kx**2)
    radius = radius / max(float(radius.max()), 1e-12)
    num_bins = max(4, min(height, width) // 2)
    bin_idx = np.clip((radius * (num_bins - 1)).astype(np.int64), 0, num_bins - 1)
    spectrum = np.bincount(bin_idx.ravel(), weights=power.ravel(), minlength=num_bins)
    counts = np.bincount(bin_idx.ravel(), minlength=num_bins)
    spectrum = spectrum / np.maximum(counts, 1)
    wavenumber = np.linspace(0.0, 1.0, num_bins)
    return wavenumber, spectrum


def spectral_band_metrics_np(
    original: np.ndarray,
    reconstructed: np.ndarray,
    *,
    edges: tuple[float, ...] = (0.0, 0.25, 0.50, 0.72, 1.01),
) -> list[dict[str, float]]:
    wave_ref, spec_ref = radial_power_spectrum_np(original)
    wave_rec, spec_rec = radial_power_spectrum_np(reconstructed)
    eps = 1e-12
    rows: list[dict[str, float]] = []
    total_ref = float(np.sum(spec_ref))
    total_rec = float(np.sum(spec_rec))
    for band_idx, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        if band_idx == len(edges) - 2:
            mask = (wave_ref >= lo) & (wave_ref <= hi)
        else:
            mask = (wave_ref >= lo) & (wave_ref < hi)
        if not np.any(mask):
            continue
        ref_band = spec_ref[mask]
        rec_band = spec_rec[mask]
        ref_power = float(np.sum(ref_band))
        rec_power = float(np.sum(rec_band))
        rows.append(
            {
                "band": float(band_idx),
                "k_min": float(lo),
                "k_max": float(hi),
                "ref_power": ref_power,
                "recon_power": rec_power,
                "power_ratio": rec_power / max(ref_power, eps),
                "ref_power_fraction": ref_power / max(total_ref, eps),
                "recon_power_fraction": rec_power / max(total_rec, eps),
                "log_l1": float(np.mean(np.abs(np.log1p(rec_band) - np.log1p(ref_band)))),
            }
        )
    return rows


def compute_variable_metric(
    original: np.ndarray,
    reconstructed: np.ndarray,
    *,
    domain: str,
    variable: str,
    sample: int,
) -> VariableMetric:
    ref = np.asarray(original, dtype=np.float64)
    pred = np.asarray(reconstructed, dtype=np.float64)
    mask = np.isfinite(ref) & np.isfinite(pred)
    if not mask.any():
        return VariableMetric(domain, variable, sample, *(float("nan"),) * 23)
    ref_v = ref[mask]
    pred_v = pred[mask]
    diff = pred_v - ref_v
    mae = float(np.mean(np.abs(diff)))
    median_ae = float(np.median(np.abs(diff)))
    p95_abs_error = float(np.percentile(np.abs(diff), 95))
    mse = float(np.mean(diff**2))
    rmse = float(math.sqrt(mse))
    bias = float(np.mean(diff))
    max_abs_error = float(np.max(np.abs(diff)))
    data_range = float(np.nanmax(ref_v) - np.nanmin(ref_v))
    ref_std = float(np.nanstd(ref_v))
    centered_diff = (pred_v - pred_v.mean()) - (ref_v - ref_v.mean())
    centered_rmse = float(math.sqrt(np.mean(centered_diff**2)))
    nmae_range = float(mae / max(data_range, 1e-12))
    nrmse_range = float(rmse / max(data_range, 1e-12))
    nrmse_std = float(rmse / max(ref_std, 1e-12))
    smape = float(100.0 * np.mean(2.0 * np.abs(diff) / np.maximum(np.abs(ref_v) + np.abs(pred_v), 1e-12)))
    mape = float(100.0 * np.mean(np.abs(diff) / np.maximum(np.abs(ref_v), 1e-12)))
    relative_l2 = float(np.linalg.norm(diff) / max(np.linalg.norm(ref_v), 1e-12))
    if ref_v.std() < 1e-12 or pred_v.std() < 1e-12:
        correlation = float("nan")
    else:
        correlation = float(np.corrcoef(ref_v, pred_v)[0, 1])
    ss_res = float(np.sum(diff**2))
    ss_tot = float(np.sum((ref_v - ref_v.mean()) ** 2))
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-12))
    explained_variance = float(1.0 - np.var(diff) / max(np.var(ref_v), 1e-12))
    nse = r2
    if np.isfinite(correlation):
        alpha = float(np.std(pred_v) / max(np.std(ref_v), 1e-12))
        beta = float(np.mean(pred_v) / max(np.mean(ref_v), 1e-12))
        kge = float(1.0 - math.sqrt((correlation - 1.0) ** 2 + (alpha - 1.0) ** 2 + (beta - 1.0) ** 2))
    else:
        kge = float("nan")
    denom = float(np.sum((np.abs(pred_v - ref_v.mean()) + np.abs(ref_v - ref_v.mean())) ** 2))
    willmott_d = float(1.0 - ss_res / max(denom, 1e-12))
    mu_x = float(ref_v.mean())
    mu_y = float(pred_v.mean())
    var_x = float(ref_v.var())
    var_y = float(pred_v.var())
    cov_xy = float(np.mean((ref_v - mu_x) * (pred_v - mu_y)))
    c1 = (0.01 * max(data_range, 1e-12)) ** 2
    c2 = (0.03 * max(data_range, 1e-12)) ** 2
    ssim_global = float(((2 * mu_x * mu_y + c1) * (2 * cov_xy + c2)) / ((mu_x**2 + mu_y**2 + c1) * (var_x + var_y + c2)))
    psnr = float(20.0 * math.log10(max(data_range, 1e-12) / max(rmse, 1e-12)))
    _, ref_spec = radial_power_spectrum_np(ref)
    _, pred_spec = radial_power_spectrum_np(pred)
    spectral_l1 = float(np.mean(np.abs(np.log1p(pred_spec) - np.log1p(ref_spec))))
    return VariableMetric(
        domain=domain,
        variable=variable,
        sample=sample,
        mae=mae,
        median_ae=median_ae,
        p95_abs_error=p95_abs_error,
        mse=mse,
        rmse=rmse,
        centered_rmse=centered_rmse,
        nmae_range=nmae_range,
        nrmse_range=nrmse_range,
        nrmse_std=nrmse_std,
        bias=bias,
        smape=smape,
        mape=mape,
        max_abs_error=max_abs_error,
        relative_l2=relative_l2,
        correlation=correlation,
        r2=r2,
        explained_variance=explained_variance,
        nse=nse,
        kge=kge,
        willmott_d=willmott_d,
        ssim_global=ssim_global,
        psnr=psnr,
        spectral_l1=spectral_l1,
    )


def metric_header() -> str:
    return (
        "domain,variable,sample,mae,median_ae,p95_abs_error,mse,rmse,centered_rmse,"
        "nmae_range,nrmse_range,nrmse_std,bias,smape,mape,max_abs_error,relative_l2,"
        "correlation,r2,explained_variance,nse,kge,willmott_d,ssim_global,psnr,spectral_l1"
    )


def metric_to_csv(metric: VariableMetric) -> str:
    values = [
        metric.domain,
        metric.variable,
        str(metric.sample),
        f"{metric.mae:.10g}",
        f"{metric.median_ae:.10g}",
        f"{metric.p95_abs_error:.10g}",
        f"{metric.mse:.10g}",
        f"{metric.rmse:.10g}",
        f"{metric.centered_rmse:.10g}",
        f"{metric.nmae_range:.10g}",
        f"{metric.nrmse_range:.10g}",
        f"{metric.nrmse_std:.10g}",
        f"{metric.bias:.10g}",
        f"{metric.smape:.10g}",
        f"{metric.mape:.10g}",
        f"{metric.max_abs_error:.10g}",
        f"{metric.relative_l2:.10g}",
        f"{metric.correlation:.10g}",
        f"{metric.r2:.10g}",
        f"{metric.explained_variance:.10g}",
        f"{metric.nse:.10g}",
        f"{metric.kge:.10g}",
        f"{metric.willmott_d:.10g}",
        f"{metric.ssim_global:.10g}",
        f"{metric.psnr:.10g}",
        f"{metric.spectral_l1:.10g}",
    ]
    return ",".join(values)


def summarize_metrics(metrics: Iterable[VariableMetric]) -> list[dict[str, float | str]]:
    grouped: dict[tuple[str, str], list[VariableMetric]] = {}
    for metric in metrics:
        grouped.setdefault((metric.domain, metric.variable), []).append(metric)
    rows: list[dict[str, float | str]] = []
    numeric_fields = (
        "mae",
        "median_ae",
        "p95_abs_error",
        "mse",
        "rmse",
        "centered_rmse",
        "nmae_range",
        "nrmse_range",
        "nrmse_std",
        "bias",
        "smape",
        "mape",
        "max_abs_error",
        "relative_l2",
        "correlation",
        "r2",
        "explained_variance",
        "nse",
        "kge",
        "willmott_d",
        "ssim_global",
        "psnr",
        "spectral_l1",
    )
    for (domain, variable), items in sorted(grouped.items()):
        row: dict[str, float | str] = {"domain": domain, "variable": variable, "count": len(items)}
        for field in numeric_fields:
            values = np.asarray([getattr(item, field) for item in items], dtype=np.float64)
            row[f"{field}_mean"] = float(np.nanmean(values))
            row[f"{field}_std"] = float(np.nanstd(values))
        rows.append(row)
    return rows


def summary_to_csv(rows: list[dict[str, float | str]]) -> str:
    if not rows:
        return "domain,variable,count\n"
    keys = list(rows[0].keys())
    lines = [",".join(keys)]
    for row in rows:
        lines.append(",".join(str(row[key]) for key in keys))
    return "\n".join(lines)


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


class DatasetView(Dataset):
    def __init__(self, dataset: Dataset, indices: list[int]) -> None:
        self.dataset = dataset
        self.indices = indices
        self.selected_indices = indices

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, idx: int):
        return self.dataset[self.indices[idx]]

    def __getattr__(self, name: str):
        return getattr(self.dataset, name)


def _strip_module_prefix(state: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
    return {key.removeprefix("module."): value for key, value in state.items()}


def _state_dict_from_checkpoint(checkpoint) -> dict[str, torch.Tensor]:
    if isinstance(checkpoint, dict):
        for key in ("model_state_dict", "state_dict", "model"):
            state = checkpoint.get(key)
            if isinstance(state, dict):
                return _strip_module_prefix(state)
        if checkpoint and all(torch.is_tensor(value) for value in checkpoint.values()):
            return _strip_module_prefix(checkpoint)
    raise ValueError("Unsupported checkpoint format. Expected a state_dict or a dict containing model_state_dict/state_dict.")


def build_model(
    *,
    img_size: list[int] | tuple[int, int] = (256, 256),
    in_chans: int = 4,
    out_chans: int = 4,
) -> LGUnet_all:
    return LGUnet_all(
        img_size=list(img_size),
        patch_size=4,
        stride=[4, 4],
        in_chans=in_chans,
        out_chans=out_chans,
        enc_depths=[2, 2, 6],
        enc_heads=[3, 6, 12],
        lg_depths=[2, 2],
        lg_heads=[6, 12],
        inchans_list=[in_chans],
        outchans_list=[out_chans],
        enc_dim=96,
        embed_dim=768,
        window_size=8,
        Weather_T=1,
        use_checkpoint=False,
        pre_norm=True,
    )


def load_model(
    checkpoint: str | Path,
    device: torch.device,
    *,
    in_chans: int = 4,
    out_chans: int = 4,
    strict: bool = True,
) -> LGUnet_all:
    ckpt = torch.load(checkpoint, map_location=device, weights_only=False)
    state = _state_dict_from_checkpoint(ckpt)
    model = build_model(in_chans=in_chans, out_chans=out_chans).to(device)
    model.load_state_dict(state, strict=strict)
    model.eval()
    return model


def _denorm(variable: str, data: torch.Tensor) -> torch.Tensor:
    base_name = variable.split("_t", 1)[0]
    stat = norm_dict[base_name]
    return data * float(stat["std"]) + float(stat["mean"])


def _write_band_metrics(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _metric_spaces(metric_space: str) -> list[str]:
    metric_space = metric_space.lower()
    if metric_space == "both":
        return ["normalized", "physical"]
    if metric_space not in {"normalized", "physical"}:
        raise ValueError("metric_space must be one of: normalized, physical, both")
    return [metric_space]


def _reconstruction_metrics_dir(output_dir: Path, metric_space: str, spaces: list[str]) -> Path:
    if len(spaces) == 1:
        return output_dir / "metrics" / "reconstruction"
    return output_dir / "metrics" / f"reconstruction_{metric_space}"


def _write_reconstruction_metrics(path: Path, metrics) -> None:
    path.mkdir(parents=True, exist_ok=True)
    summary = summarize_metrics(metrics)
    (path / "per_sample_variable_metrics.csv").write_text(
        metric_header() + "\n" + "\n".join(metric_to_csv(item) for item in metrics),
        encoding="utf-8",
    )
    (path / "summary_metrics.csv").write_text(summary_to_csv(summary), encoding="utf-8")
    (path / "summary_metrics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def _masked_array(data: torch.Tensor, mask: torch.Tensor | None) -> tuple[np.ndarray, np.ndarray]:
    arr = data.detach().cpu().numpy()
    if mask is None:
        return arr, arr
    valid = mask.detach().cpu().numpy().astype(bool)
    return np.where(valid, arr, np.nan), np.where(valid, arr, 0.0)


def resolve_analysis_output_dir(checkpoint: str | Path, output_dir: str | Path | None = None) -> Path:
    if output_dir is not None:
        return Path(output_dir)
    checkpoint_path = Path(checkpoint)
    run_dir = checkpoint_path.parent if checkpoint_path.suffix else checkpoint_path
    return run_dir / "analysis_result"


def _load_config(path: str | Path | None) -> dict:
    if path is None:
        return {}
    config_path = Path(path)
    with open(config_path, "r", encoding="utf-8") as f:
        text = f.read()
    if config_path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError("YAML analysis configs require PyYAML; use JSON or install PyYAML.") from exc
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    return data or {}


def _as_int_list(values) -> list[int] | None:
    if values is None:
        return None
    return [int(item) for item in values]


def _variable_names(variables: list[str], length: int) -> list[str]:
    if length <= 1:
        return list(variables)
    return [f"{variable}_t{step}" for step in range(length) for variable in variables]


def build_analysis_dataset(
    data_dir: str | Path,
    *,
    input_len: int = 1,
    label_len: int = 1,
    input_vars: list[str] | None = None,
    label_vars: list[str] | None = None,
    split: str = "val",
    start_month: int = 1,
    end_month: int = 12,
    base_year: int = 2022,
) -> Dataset:
    input_vars = input_vars or list(DEFAULT_VARIABLES)
    label_vars = label_vars or list(DEFAULT_VARIABLES)
    dataset = WRFSlidingDataset(
        nc_dir=str(data_dir),
        input_len=input_len,
        label_len=label_len,
        input_vars=input_vars,
        label_vars=label_vars,
        start_month=start_month,
        end_month=end_month,
    )
    dataset.analysis_variables = _variable_names(label_vars, label_len)
    split = split.lower()
    if split == "all":
        dataset.selected_indices = list(range(len(dataset)))
        return dataset
    quarter_ids = dataset.get_sample_quarters(base_year=base_year)
    if split == "train":
        indices = [idx for idx, quarter in enumerate(quarter_ids) if 1 <= quarter <= 10]
    elif split == "val":
        indices = [idx for idx, quarter in enumerate(quarter_ids) if quarter == 11]
    elif split == "test":
        indices = [idx for idx, quarter in enumerate(quarter_ids) if quarter == 12]
    else:
        raise ValueError("split must be one of: train, val, test, all")
    view = DatasetView(dataset, indices)
    view.analysis_variables = dataset.analysis_variables
    return view


def _select_analysis_samples(
    dataset,
    *,
    sample_start: int = 0,
    sample_stop: int | None = None,
    sample_indices: list[int] | None = None,
):
    if sample_indices is not None:
        indices = [idx for idx in sample_indices if 0 <= idx < len(dataset)]
    else:
        start = max(int(sample_start), 0)
        stop = len(dataset) if sample_stop is None else min(int(sample_stop), len(dataset))
        indices = list(range(start, max(start, stop)))
    view = DatasetView(dataset, indices)
    view.analysis_variables = getattr(dataset, "analysis_variables", list(DEFAULT_VARIABLES))
    source = list(getattr(dataset, "selected_indices", range(len(dataset))))
    view.selected_indices = [int(source[idx]) if idx < len(source) else int(idx) for idx in indices]
    return view


def _unpack_model_output(output):
    if isinstance(output, dict):
        for key in ("recon", "reconstruction", "output", "pred", "prediction"):
            if key in output:
                return output[key]
        raise ValueError(f"Model returned a dict without a known reconstruction key: {list(output.keys())}")
    if isinstance(output, (tuple, list)):
        return output[0]
    return output


@torch.no_grad()
def run_analysis(
    checkpoint: str | Path = DEFAULT_MODEL_PATH,
    data_dir: str | Path = DEFAULT_DATA_DIR,
    output_dir: str | Path | None = None,
    num_samples: int | None = 4,
    max_variables: int | None = None,
    plot_space: str = "normalized",
    strict_load: bool = True,
    split: str = "val",
    sample_start: int = 0,
    sample_stop: int | None = None,
    sample_indices: list[int] | None = None,
    input_len: int = 1,
    label_len: int = 1,
    input_vars: list[str] | None = None,
    label_vars: list[str] | None = None,
    start_month: int = 1,
    end_month: int = 12,
    base_year: int = 2022,
    cmap: str = "bwr",
    spectrum_colors: tuple[str, str] = ("#1f77b4", "#d62728"),
    make_plots: bool = True,
    metric_space: str = "physical",
) -> Path:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    input_vars = input_vars or list(DEFAULT_VARIABLES)
    label_vars = label_vars or list(DEFAULT_VARIABLES)
    dataset = build_analysis_dataset(
        data_dir,
        input_len=input_len,
        label_len=label_len,
        input_vars=input_vars,
        label_vars=label_vars,
        split=split,
        start_month=start_month,
        end_month=end_month,
        base_year=base_year,
    )
    dataset = _select_analysis_samples(
        dataset,
        sample_start=sample_start,
        sample_stop=sample_stop,
        sample_indices=sample_indices,
    )
    in_chans = len(input_vars) * input_len
    out_chans = len(label_vars) * label_len
    model = load_model(checkpoint, device, in_chans=in_chans, out_chans=out_chans, strict=strict_load)
    output_dir = resolve_analysis_output_dir(checkpoint, output_dir)
    reconstruction_dir = output_dir / "plots" / "reconstruction"
    spectrum_dir = output_dir / "plots" / "spectra"
    spectral_metrics_dir = output_dir / "metrics" / "spectral_bands"
    if make_plots:
        reconstruction_dir.mkdir(parents=True, exist_ok=True)
        spectrum_dir.mkdir(parents=True, exist_ok=True)
    spectral_metrics_dir.mkdir(parents=True, exist_ok=True)

    spaces = _metric_spaces(metric_space)
    metrics_by_space = {space: [] for space in spaces}
    band_rows: list[dict[str, float | int | str]] = []
    plot_space = plot_space.lower()
    if plot_space not in {"normalized", "physical", "both"}:
        raise ValueError("plot_space must be one of: normalized, physical, both")
    sample_limit = len(dataset) if num_samples is None or int(num_samples) <= 0 else min(int(num_samples), len(dataset))
    source_indices = list(getattr(dataset, "selected_indices", range(len(dataset))))
    variables = list(getattr(dataset, "analysis_variables", _variable_names(label_vars, label_len)))
    for idx in range(sample_limit):
        source_idx = int(source_indices[idx]) if idx < len(source_indices) else idx
        sample = dataset[idx]
        if len(sample) == 3:
            inputs, targets, mask = sample
        else:
            inputs, targets = sample
            mask = torch.ones_like(targets)
        output = _unpack_model_output(model(inputs.unsqueeze(0).to(device)))
        recon_tensor = output.squeeze(0).detach().cpu()
        target_tensor = targets.detach().cpu()
        mask_tensor = mask.detach().cpu() if mask is not None else torch.ones_like(target_tensor)

        limit = len(variables) if max_variables is None else min(max_variables, len(variables))
        for channel, variable in enumerate(variables[:limit]):
            channel_mask = mask_tensor[channel] if mask_tensor is not None else None
            original_norm, original_norm_fft = _masked_array(target_tensor[channel], channel_mask)
            reconstructed_norm, reconstructed_norm_fft = _masked_array(recon_tensor[channel], channel_mask)
            original_phys, _ = _masked_array(_denorm(variable, target_tensor[channel]), channel_mask)
            reconstructed_phys, _ = _masked_array(_denorm(variable, recon_tensor[channel]), channel_mask)
            if "normalized" in metrics_by_space:
                metrics_by_space["normalized"].append(
                    compute_variable_metric(
                        original_norm,
                        reconstructed_norm,
                        domain="wrf",
                        variable=variable,
                        sample=idx,
                    )
                )
            if "physical" in metrics_by_space:
                metrics_by_space["physical"].append(
                    compute_variable_metric(
                        original_phys,
                        reconstructed_phys,
                        domain="wrf",
                        variable=variable,
                        sample=idx,
                    )
                )
            stem = f"sample_{idx:04d}_wrf_{variable}"
            if make_plots:
                plot_pairs = []
                if plot_space in {"normalized", "both"}:
                    plot_pairs.append(("normalized", original_norm, reconstructed_norm))
                if plot_space in {"physical", "both"}:
                    plot_pairs.append(("physical", original_phys, reconstructed_phys))
                for suffix, original, reconstructed in plot_pairs:
                    name = stem if plot_space != "both" else f"{stem}_{suffix}"
                    title_suffix = "normalized" if suffix == "normalized" else "physical"
                    plot_reconstruction_triplet(
                        original,
                        reconstructed,
                        reconstruction_dir / f"{name}_triptych.png",
                        f"WRF {variable} frame {source_idx} ({title_suffix})",
                        cmap=cmap,
                    )
                    if suffix == "normalized":
                        spectrum_original = original_norm_fft
                        spectrum_reconstructed = reconstructed_norm_fft
                    else:
                        spectrum_original = np.nan_to_num(original, nan=0.0)
                        spectrum_reconstructed = np.nan_to_num(reconstructed, nan=0.0)
                    wave_ref, spec_ref = radial_power_spectrum_np(spectrum_original)
                    _, spec_rec = radial_power_spectrum_np(spectrum_reconstructed)
                    plot_spectrum_pair(
                        wave_ref,
                        spec_ref,
                        spec_rec,
                        spectrum_dir / f"{name}_spectrum.png",
                        f"WRF {variable} radial spectrum frame {source_idx} ({title_suffix})",
                        colors=spectrum_colors,
                    )
            for row in spectral_band_metrics_np(original_norm_fft, reconstructed_norm_fft):
                band_rows.append(
                    {
                        "domain": "wrf",
                        "variable": variable,
                        "sample": idx,
                        **row,
                    }
                )

    for space, metrics in metrics_by_space.items():
        _write_reconstruction_metrics(_reconstruction_metrics_dir(output_dir, space, spaces), metrics)
    _write_band_metrics(spectral_metrics_dir / "spectral_band_metrics.csv", band_rows)
    (spectral_metrics_dir / "spectral_band_metrics.json").write_text(json.dumps(band_rows, indent=2), encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "analysis.json"))
    parser.add_argument("--checkpoint", default=str(DEFAULT_MODEL_PATH))
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR))
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--num-samples", type=int, default=None)
    parser.add_argument("--all-samples", action="store_true")
    parser.add_argument("--max-variables", type=int, default=None)
    parser.add_argument("--plot-space", choices=["normalized", "physical", "both"], default="normalized")
    parser.add_argument("--metric-space", choices=["normalized", "physical", "both"], default=None)
    parser.add_argument("--split", choices=["train", "val", "test", "all"], default=None)
    parser.add_argument("--sample-start", type=int, default=None)
    parser.add_argument("--sample-stop", type=int, default=None)
    parser.add_argument("--sample-indices", nargs="+", type=int, default=None)
    parser.add_argument("--input-len", type=int, default=None)
    parser.add_argument("--label-len", type=int, default=None)
    parser.add_argument("--input-vars", nargs="+", default=None)
    parser.add_argument("--label-vars", nargs="+", default=None)
    parser.add_argument("--start-month", type=int, default=None)
    parser.add_argument("--end-month", type=int, default=None)
    parser.add_argument("--base-year", type=int, default=None)
    parser.add_argument("--cmap", default=None)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--non-strict-load", action="store_true")
    args = parser.parse_args()

    config = _load_config(args.config) if args.config and Path(args.config).exists() else {}
    spectrum_colors = config.get("spectrum_colors", ["#1f77b4", "#d62728"])
    config_all_samples = bool(config.get("all_samples", False))
    if args.all_samples or config_all_samples:
        num_samples = None
    else:
        num_samples = args.num_samples if args.num_samples is not None else config.get("num_samples", 4)
    run_analysis(
        args.checkpoint if args.checkpoint != str(DEFAULT_MODEL_PATH) else config.get("checkpoint", args.checkpoint),
        args.data_dir if args.data_dir != str(DEFAULT_DATA_DIR) else config.get("data_dir", args.data_dir),
        args.output_dir if args.output_dir is not None else config.get("output_dir"),
        num_samples,
        args.max_variables if args.max_variables is not None else config.get("max_variables"),
        plot_space=args.plot_space if args.plot_space != "normalized" else config.get("plot_space", args.plot_space),
        strict_load=not args.non_strict_load,
        split=args.split or config.get("split", "val"),
        sample_start=args.sample_start if args.sample_start is not None else int(config.get("sample_start", 0)),
        sample_stop=args.sample_stop if args.sample_stop is not None else config.get("sample_stop"),
        sample_indices=args.sample_indices if args.sample_indices is not None else _as_int_list(config.get("sample_indices")),
        input_len=args.input_len if args.input_len is not None else int(config.get("input_len", 1)),
        label_len=args.label_len if args.label_len is not None else int(config.get("label_len", 1)),
        input_vars=args.input_vars if args.input_vars is not None else config.get("input_vars"),
        label_vars=args.label_vars if args.label_vars is not None else config.get("label_vars"),
        start_month=args.start_month if args.start_month is not None else int(config.get("start_month", 1)),
        end_month=args.end_month if args.end_month is not None else int(config.get("end_month", 12)),
        base_year=args.base_year if args.base_year is not None else int(config.get("base_year", 2022)),
        cmap=args.cmap or config.get("cmap", "bwr"),
        spectrum_colors=(str(spectrum_colors[0]), str(spectrum_colors[1])),
        make_plots=not (args.no_plot or bool(config.get("no_plot", False))),
        metric_space=args.metric_space or config.get("metric_space", "physical"),
    )


if __name__ == "__main__":
    main()
