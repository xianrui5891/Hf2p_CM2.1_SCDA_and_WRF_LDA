from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import torch


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
    """Compare radial spectrum power by normalized wavenumber bands."""
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
    compute_spectral: bool = True,
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
    if compute_spectral:
        _, ref_spec = radial_power_spectrum_np(ref)
        _, pred_spec = radial_power_spectrum_np(pred)
        spectral_l1 = float(np.mean(np.abs(np.log1p(pred_spec) - np.log1p(ref_spec))))
    else:
        spectral_l1 = float("nan")
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
            if np.isfinite(values).any():
                row[f"{field}_mean"] = float(np.nanmean(values))
                row[f"{field}_std"] = float(np.nanstd(values))
            else:
                row[f"{field}_mean"] = float("nan")
                row[f"{field}_std"] = float("nan")
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
