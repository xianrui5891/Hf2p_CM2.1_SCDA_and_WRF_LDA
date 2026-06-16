from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from analysis.metrics import (
    VariableMetric,
    compute_variable_metric,
    metric_header,
    metric_to_csv,
    radial_power_spectrum_np,
    spectral_band_metrics_np,
    summarize_metrics,
    summary_to_csv,
)
from model import CoupledAE
from utils.data import CM2DatasetView, build_dataset
from utils.paths import DEFAULT_AE_MODEL_PATH, DEFAULT_METADATA_PATH, DEFAULT_PROCESSED_PT
from utils.plotting import plot_reconstruction_triplet, plot_spectrum_pair


def _infer_fusion_mode(cfg: dict, state: dict[str, torch.Tensor]) -> str:
    if cfg.get("fusion_mode"):
        return str(cfg["fusion_mode"])
    coupling_keys = [key for key in state if ".pre_down_coupling." in key or ".post_down_coupling." in key]
    if any(".net." in key for key in coupling_keys):
        return "cnn"
    return "cross_attention"


def load_model(checkpoint: str | Path, device: torch.device, strict: bool = True) -> CoupledAE:
    ckpt = torch.load(checkpoint, map_location=device, weights_only=False)
    cfg = ckpt.get("model_config", {})
    state = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt))
    state = {key.removeprefix("module."): value for key, value in state.items()}
    width = int(cfg.get("width", 64))
    model = CoupledAE(
        latent_dim=cfg.get("latent_dim"),
        latent_channels=cfg.get("latent_channels"),
        in_atm_channels=cfg.get("in_atm_channels", 4),
        in_ocn_channels=cfg.get("in_ocn_channels", 4),
        width=width,
        decoder_width=cfg.get("decoder_width", width * 2),
        decoder_middle_depth=cfg.get("decoder_middle_depth", 2),
        decoder_stage_depth=cfg.get("decoder_stage_depth", 1),
        atm_decoder_stage_depth=cfg.get("atm_decoder_stage_depth"),
        ocn_decoder_stage_depth=cfg.get("ocn_decoder_stage_depth"),
        cross_depth=cfg.get("cross_depth", 4),
        heads=cfg.get("heads", 4),
        attention_backend=cfg.get("attention_backend", "auto"),
        fusion_mode=_infer_fusion_mode(cfg, state),
        cnn_fusion_depth=cfg.get("cnn_fusion_depth", 2),
    ).to(device)
    model.load_state_dict(state, strict=strict)
    model.eval()
    return model


def _denorm(dataset, domain: str, variable: str, data: torch.Tensor) -> torch.Tensor:
    stat = dataset.normalization_stats[f"{domain}::{variable}"]
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


def _masked_array(data: torch.Tensor, mask: torch.Tensor | None) -> tuple:
    arr = data.numpy()
    if mask is None:
        return arr, arr
    valid = mask.numpy().astype(bool)
    return np.where(valid, arr, np.nan), np.where(valid, arr, 0.0)


def _nan_stat(values: np.ndarray, fn, default: float = float("nan")) -> np.ndarray:
    with np.errstate(all="ignore"):
        result = fn(values, axis=1)
    result = np.asarray(result, dtype=np.float64)
    result[~np.isfinite(result)] = default
    return result


def _compute_variable_metrics_batch(
    original: torch.Tensor,
    reconstructed: torch.Tensor,
    mask: torch.Tensor | None,
    *,
    domain: str,
    variable: str,
    sample_offset: int,
    normalized: bool,
) -> list[VariableMetric]:
    ref = original.detach().cpu().numpy().astype(np.float64, copy=False)
    pred = reconstructed.detach().cpu().numpy().astype(np.float64, copy=False)
    if mask is None:
        valid = np.isfinite(ref) & np.isfinite(pred)
    else:
        valid = mask.detach().cpu().numpy().astype(bool) & np.isfinite(ref) & np.isfinite(pred)
    batch = int(ref.shape[0])
    ref_flat = ref.reshape(batch, -1)
    pred_flat = pred.reshape(batch, -1)
    valid_flat = valid.reshape(batch, -1)
    count = valid_flat.sum(axis=1).astype(np.float64)
    safe_count = np.maximum(count, 1.0)
    ref_nan = np.where(valid_flat, ref_flat, np.nan)
    pred_nan = np.where(valid_flat, pred_flat, np.nan)
    diff_nan = pred_nan - ref_nan
    abs_diff = np.abs(diff_nan)
    diff_zero = np.nan_to_num(diff_nan, nan=0.0)
    ref_zero = np.nan_to_num(ref_nan, nan=0.0)
    pred_zero = np.nan_to_num(pred_nan, nan=0.0)
    abs_diff_zero = np.nan_to_num(abs_diff, nan=0.0)

    mae = abs_diff_zero.sum(axis=1) / safe_count
    median_ae = _nan_stat(abs_diff, np.nanmedian)
    p95_abs_error = _nan_stat(abs_diff, lambda arr, axis: np.nanpercentile(arr, 95, axis=axis))
    mse = (diff_zero**2).sum(axis=1) / safe_count
    rmse = np.sqrt(mse)
    bias = diff_zero.sum(axis=1) / safe_count
    max_abs_error = _nan_stat(abs_diff, np.nanmax)
    ref_min = _nan_stat(ref_nan, np.nanmin)
    ref_max = _nan_stat(ref_nan, np.nanmax)
    data_range = ref_max - ref_min
    ref_mean = ref_zero.sum(axis=1) / safe_count
    pred_mean = pred_zero.sum(axis=1) / safe_count
    ref_center = np.where(valid_flat, ref_flat - ref_mean[:, None], 0.0)
    pred_center = np.where(valid_flat, pred_flat - pred_mean[:, None], 0.0)
    ref_var = (ref_center**2).sum(axis=1) / safe_count
    pred_var = (pred_center**2).sum(axis=1) / safe_count
    ref_std = np.sqrt(ref_var)
    pred_std = np.sqrt(pred_var)
    centered_rmse = np.sqrt(((pred_center - ref_center) ** 2).sum(axis=1) / safe_count)
    nmae_range = mae / np.maximum(data_range, 1e-12)
    nrmse_range = rmse / np.maximum(data_range, 1e-12)
    nrmse_std = rmse / np.maximum(ref_std, 1e-12)
    smape = 100.0 * np.nansum(2.0 * abs_diff / np.maximum(np.abs(ref_nan) + np.abs(pred_nan), 1e-12), axis=1) / safe_count
    mape = 100.0 * np.nansum(abs_diff / np.maximum(np.abs(ref_nan), 1e-12), axis=1) / safe_count
    relative_l2 = np.sqrt((diff_zero**2).sum(axis=1)) / np.maximum(np.sqrt((ref_zero**2).sum(axis=1)), 1e-12)
    cov = (ref_center * pred_center).sum(axis=1) / safe_count
    correlation = cov / np.maximum(ref_std * pred_std, 1e-12)
    correlation[(ref_std < 1e-12) | (pred_std < 1e-12)] = np.nan
    ss_res = (diff_zero**2).sum(axis=1)
    ss_tot = (ref_center**2).sum(axis=1)
    r2 = 1.0 - ss_res / np.maximum(ss_tot, 1e-12)
    explained_variance = 1.0 - np.maximum(((diff_zero - bias[:, None]) ** 2 * valid_flat).sum(axis=1) / safe_count, 0.0) / np.maximum(ref_var, 1e-12)
    nse = r2
    alpha = pred_std / np.maximum(ref_std, 1e-12)
    beta = pred_mean / np.maximum(ref_mean, 1e-12)
    kge = 1.0 - np.sqrt((correlation - 1.0) ** 2 + (alpha - 1.0) ** 2 + (beta - 1.0) ** 2)
    kge[~np.isfinite(correlation)] = np.nan
    denom = ((np.abs(pred_flat - ref_mean[:, None]) + np.abs(ref_flat - ref_mean[:, None])) ** 2)
    denom = np.where(valid_flat, denom, 0.0).sum(axis=1)
    willmott_d = 1.0 - ss_res / np.maximum(denom, 1e-12)
    c1 = (0.01 * np.maximum(data_range, 1e-12)) ** 2
    c2 = (0.03 * np.maximum(data_range, 1e-12)) ** 2
    ssim_global = ((2 * ref_mean * pred_mean + c1) * (2 * cov + c2)) / (
        (ref_mean**2 + pred_mean**2 + c1) * (ref_var + pred_var + c2)
    )
    psnr = 20.0 * np.log10(np.maximum(data_range, 1e-12) / np.maximum(rmse, 1e-12))
    spectral_l1 = np.full(batch, np.nan, dtype=np.float64)

    metrics: list[VariableMetric] = []
    for item in range(batch):
        if count[item] < 1:
            metrics.append(VariableMetric(domain, variable, sample_offset + item, *(float("nan"),) * 23))
            continue
        metrics.append(
            VariableMetric(
                domain=domain,
                variable=variable,
                sample=sample_offset + item,
                mae=float(mae[item]),
                median_ae=float(median_ae[item]),
                p95_abs_error=float(p95_abs_error[item]),
                mse=float(mse[item]),
                rmse=float(rmse[item]),
                centered_rmse=float(centered_rmse[item]),
                nmae_range=float(nmae_range[item]),
                nrmse_range=float(nrmse_range[item]),
                nrmse_std=float(nrmse_std[item]),
                bias=float(bias[item]),
                smape=float(smape[item]),
                mape=float(mape[item]),
                max_abs_error=float(max_abs_error[item]),
                relative_l2=float(relative_l2[item]),
                correlation=float(correlation[item]),
                r2=float(r2[item]),
                explained_variance=float(explained_variance[item]),
                nse=float(nse[item]),
                kge=float(kge[item]),
                willmott_d=float(willmott_d[item]),
                ssim_global=float(ssim_global[item]),
                psnr=float(psnr[item]),
                spectral_l1=float(spectral_l1[item]),
            )
        )
    return metrics


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
        except ImportError as exc:  # pragma: no cover - only used when yaml configs are requested.
            raise RuntimeError("YAML analysis configs require PyYAML; use JSON or install PyYAML.") from exc
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    return data or {}


def _as_int_list(values) -> list[int] | None:
    if values is None:
        return None
    return [int(item) for item in values]


def _progress_iter(iterable, *, total: int, enabled: bool, desc: str):
    if not enabled:
        return iterable
    try:
        from tqdm.auto import tqdm

        return tqdm(iterable, total=total, desc=desc)
    except ImportError:
        step = max(1, total // 20)

        def _printer():
            for idx, item in enumerate(iterable, 1):
                if idx == 1 or idx == total or idx % step == 0:
                    print(f"{desc}: {idx}/{total} batches", flush=True)
                yield item

        return _printer()


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
    return CM2DatasetView(dataset, indices)


@torch.no_grad()
def run_analysis(
    checkpoint: str | Path = DEFAULT_AE_MODEL_PATH,
    data_path: str | Path = DEFAULT_PROCESSED_PT,
    metadata_path: str | Path = DEFAULT_METADATA_PATH,
    output_dir: str | Path | None = None,
    num_samples: int | None = 4,
    max_variables: int | None = None,
    plot_space: str = "normalized",
    strict_load: bool = True,
    split: str = "val",
    val_fraction: float = 0.1,
    sample_start: int = 0,
    sample_stop: int | None = None,
    sample_indices: list[int] | None = None,
    cmap: str = "bwr",
    spectrum_colors: tuple[str, str] = ("#1f77b4", "#d62728"),
    make_plots: bool = True,
    metric_space: str = "physical",
    batch_size: int = 1,
    num_workers: int = 0,
    show_progress: bool = True,
    include_spectral: bool = True,
) -> Path:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dataset = build_dataset(data_path, metadata_path=metadata_path, split=split, val_fraction=val_fraction)
    dataset = _select_analysis_samples(
        dataset,
        sample_start=sample_start,
        sample_stop=sample_stop,
        sample_indices=sample_indices,
    )
    model = load_model(checkpoint, device, strict=strict_load)
    output_dir = resolve_analysis_output_dir(checkpoint, output_dir)
    reconstruction_dir = output_dir / "plots" / "reconstruction"
    spectrum_dir = output_dir / "plots" / "spectra"
    spectral_metrics_dir = output_dir / "metrics" / "spectral_bands"
    if make_plots:
        reconstruction_dir.mkdir(parents=True, exist_ok=True)
    if make_plots and include_spectral:
        spectrum_dir.mkdir(parents=True, exist_ok=True)
    if include_spectral:
        spectral_metrics_dir.mkdir(parents=True, exist_ok=True)

    spaces = _metric_spaces(metric_space)
    metrics_by_space = {space: [] for space in spaces}
    band_rows: list[dict[str, float | int | str]] = []
    plot_space = plot_space.lower()
    if plot_space not in {"normalized", "physical", "both"}:
        raise ValueError("plot_space must be one of: normalized, physical, both")
    sample_limit = len(dataset) if num_samples is None or int(num_samples) <= 0 else min(int(num_samples), len(dataset))
    source_indices = list(getattr(dataset, "selected_indices", range(len(dataset))))
    source_indices = source_indices[:sample_limit]
    analysis_dataset = _select_analysis_samples(dataset, sample_indices=list(range(sample_limit)))
    batch_size = max(1, int(batch_size))
    loader = DataLoader(
        analysis_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=max(0, int(num_workers)),
        pin_memory=device.type == "cuda",
    )
    progress = _progress_iter(loader, total=len(loader), enabled=show_progress, desc=f"Analyzing {split}")
    processed = 0
    for batch in progress:
        if len(batch) == 4:
            atm_batch, ocn_batch, atm_mask_batch, ocn_mask_batch = batch
        else:
            atm_batch, ocn_batch = batch
            atm_mask_batch = torch.ones_like(atm_batch, dtype=torch.bool)
            ocn_mask_batch = torch.ones_like(ocn_batch, dtype=torch.bool)
        out = model(atm_batch.to(device, non_blocking=True), ocn_batch.to(device, non_blocking=True))
        atm_recon_batch = out["atm_recon"].detach().cpu()
        ocn_recon_batch = out["ocn_recon"].detach().cpu()
        batch_count = int(atm_batch.shape[0])
        if not make_plots and not include_spectral:
            batch_domains = [
                ("atm", dataset.atm_vars, atm_batch, atm_recon_batch, atm_mask_batch),
                ("ocn", dataset.ocn_vars, ocn_batch, ocn_recon_batch, ocn_mask_batch),
            ]
            for domain, variables, original_batch, recon_batch, mask_batch in batch_domains:
                limit = len(variables) if max_variables is None else min(max_variables, len(variables))
                for channel, variable in enumerate(variables[:limit]):
                    channel_mask_batch = mask_batch[:, channel] if mask_batch is not None else None
                    if "normalized" in metrics_by_space:
                        metrics_by_space["normalized"].extend(
                            _compute_variable_metrics_batch(
                                original_batch[:, channel],
                                recon_batch[:, channel],
                                channel_mask_batch,
                                domain=domain,
                                variable=variable,
                                sample_offset=processed,
                                normalized=True,
                            )
                        )
                    if "physical" in metrics_by_space:
                        metrics_by_space["physical"].extend(
                            _compute_variable_metrics_batch(
                                _denorm(dataset, domain, variable, original_batch[:, channel]),
                                _denorm(dataset, domain, variable, recon_batch[:, channel]),
                                channel_mask_batch,
                                domain=domain,
                                variable=variable,
                                sample_offset=processed,
                                normalized=False,
                            )
                        )
            processed += batch_count
            continue
        for batch_idx in range(batch_count):
            idx = processed + batch_idx
            source_idx = int(source_indices[idx]) if idx < len(source_indices) else idx
            domains = [
                ("atm", dataset.atm_vars, atm_batch[batch_idx], atm_recon_batch[batch_idx], atm_mask_batch[batch_idx]),
                ("ocn", dataset.ocn_vars, ocn_batch[batch_idx], ocn_recon_batch[batch_idx], ocn_mask_batch[batch_idx]),
            ]
            for domain, variables, original_tensor, recon_tensor, mask_tensor in domains:
                limit = len(variables) if max_variables is None else min(max_variables, len(variables))
                for channel, variable in enumerate(variables[:limit]):
                    channel_mask = mask_tensor[channel] if mask_tensor is not None else None
                    original_norm, original_norm_fft = _masked_array(original_tensor[channel], channel_mask)
                    reconstructed_norm, reconstructed_norm_fft = _masked_array(recon_tensor[channel], channel_mask)
                    original_phys, _ = _masked_array(_denorm(dataset, domain, variable, original_tensor[channel]), channel_mask)
                    reconstructed_phys, _ = _masked_array(_denorm(dataset, domain, variable, recon_tensor[channel]), channel_mask)
                    if "normalized" in metrics_by_space:
                        metrics_by_space["normalized"].append(
                            compute_variable_metric(
                                original_norm,
                                reconstructed_norm,
                                domain=domain,
                                variable=variable,
                                sample=idx,
                                compute_spectral=include_spectral,
                            )
                        )
                    if "physical" in metrics_by_space:
                        metrics_by_space["physical"].append(
                            compute_variable_metric(
                                original_phys,
                                reconstructed_phys,
                                domain=domain,
                                variable=variable,
                                sample=idx,
                                compute_spectral=include_spectral,
                            )
                        )
                    stem = f"sample_{idx:04d}_{domain}_{variable}"
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
                                f"{domain.upper()} {variable} frame {source_idx} ({title_suffix})",
                                cmap=cmap,
                            )
                            if include_spectral:
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
                                    f"{domain.upper()} {variable} radial spectrum frame {source_idx} ({title_suffix})",
                                    colors=spectrum_colors,
                                )
                    if include_spectral:
                        for row in spectral_band_metrics_np(original_norm_fft, reconstructed_norm_fft):
                            band_rows.append(
                                {
                                    "domain": domain,
                                    "variable": variable,
                                    "sample": idx,
                                    **row,
                                }
                            )
        processed += batch_count

    for space, metrics in metrics_by_space.items():
        _write_reconstruction_metrics(_reconstruction_metrics_dir(output_dir, space, spaces), metrics)
    if include_spectral:
        _write_band_metrics(spectral_metrics_dir / "spectral_band_metrics.csv", band_rows)
        (spectral_metrics_dir / "spectral_band_metrics.json").write_text(json.dumps(band_rows, indent=2), encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "analysis.json"))
    parser.add_argument("--checkpoint", default=str(DEFAULT_AE_MODEL_PATH))
    parser.add_argument("--data", default=str(DEFAULT_PROCESSED_PT))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA_PATH))
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--num-samples", type=int, default=None)
    parser.add_argument("--all-samples", action="store_true")
    parser.add_argument("--max-variables", type=int, default=None)
    parser.add_argument("--plot-space", choices=["normalized", "physical", "both"], default="normalized")
    parser.add_argument("--metric-space", choices=["normalized", "physical", "both"], default=None)
    parser.add_argument("--split", choices=["train", "val", "all"], default=None)
    parser.add_argument("--val-only", action="store_true")
    parser.add_argument("--val-fraction", type=float, default=None)
    parser.add_argument("--sample-start", type=int, default=None)
    parser.add_argument("--sample-stop", type=int, default=None)
    parser.add_argument("--sample-indices", nargs="+", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--num-workers", type=int, default=None)
    parser.add_argument("--cmap", default=None)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--no-spectral", action="store_true")
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--non-strict-load", action="store_true")
    args = parser.parse_args()
    config = _load_config(args.config) if args.config and Path(args.config).exists() else {}
    spectrum_colors = config.get("spectrum_colors", ["#1f77b4", "#d62728"])
    config_all_samples = bool(config.get("all_samples", False))
    if args.all_samples or config_all_samples:
        num_samples = None
    else:
        num_samples = args.num_samples if args.num_samples is not None else config.get("num_samples", 4)
    split = "val" if args.val_only else (args.split or config.get("split", "val"))
    run_analysis(
        args.checkpoint if args.checkpoint != str(DEFAULT_AE_MODEL_PATH) else config.get("checkpoint", args.checkpoint),
        args.data if args.data != str(DEFAULT_PROCESSED_PT) else config.get("data_path", args.data),
        args.metadata if args.metadata != str(DEFAULT_METADATA_PATH) else config.get("metadata_path", args.metadata),
        args.output_dir if args.output_dir is not None else config.get("output_dir"),
        num_samples,
        args.max_variables if args.max_variables is not None else config.get("max_variables"),
        plot_space=args.plot_space if args.plot_space != "normalized" else config.get("plot_space", args.plot_space),
        strict_load=not args.non_strict_load,
        split=split,
        val_fraction=args.val_fraction if args.val_fraction is not None else float(config.get("val_fraction", 0.1)),
        sample_start=args.sample_start if args.sample_start is not None else int(config.get("sample_start", 0)),
        sample_stop=args.sample_stop if args.sample_stop is not None else config.get("sample_stop"),
        sample_indices=args.sample_indices if args.sample_indices is not None else _as_int_list(config.get("sample_indices")),
        cmap=args.cmap or config.get("cmap", "bwr"),
        spectrum_colors=(str(spectrum_colors[0]), str(spectrum_colors[1])),
        make_plots=not (args.no_plot or bool(config.get("no_plot", False))),
        metric_space=args.metric_space or config.get("metric_space", "physical"),
        batch_size=args.batch_size if args.batch_size is not None else int(config.get("batch_size", 1)),
        num_workers=args.num_workers if args.num_workers is not None else int(config.get("num_workers", 0)),
        show_progress=not (args.no_progress or bool(config.get("no_progress", False))),
        include_spectral=not (args.no_spectral or bool(config.get("no_spectral", False))),
    )


if __name__ == "__main__":
    main()
