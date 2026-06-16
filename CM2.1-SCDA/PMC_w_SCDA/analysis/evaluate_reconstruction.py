from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from analysis.analyze_ae import load_model, resolve_analysis_output_dir
from analysis.metrics import compute_variable_metric, metric_header, metric_to_csv, summarize_metrics, summary_to_csv
from utils.data import build_dataset
from utils.paths import DEFAULT_AE_MODEL_PATH, DEFAULT_METADATA_PATH, DEFAULT_PROCESSED_PT


def _denorm(dataset, domain: str, variable: str, data: torch.Tensor) -> torch.Tensor:
    stat = dataset.normalization_stats[f"{domain}::{variable}"]
    return data * float(stat["std"]) + float(stat["mean"])


def _masked_numpy(data: torch.Tensor, mask: torch.Tensor | None) -> np.ndarray:
    arr = data.numpy()
    if mask is None:
        return arr
    return np.where(mask.numpy().astype(bool), arr, np.nan)


@torch.no_grad()
def evaluate(
    checkpoint: str | Path = DEFAULT_AE_MODEL_PATH,
    data_path: str | Path = DEFAULT_PROCESSED_PT,
    metadata_path: str | Path = DEFAULT_METADATA_PATH,
    output_dir: str | Path | None = None,
    batch_size: int = 2,
    max_samples: int | None = None,
    split: str = "val",
    val_fraction: float = 0.1,
) -> Path:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dataset = build_dataset(
        data_path,
        metadata_path=metadata_path,
        split=split,
        val_fraction=val_fraction,
        max_samples=max_samples,
    )
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)
    model = load_model(checkpoint, device)
    output_dir = resolve_analysis_output_dir(checkpoint, output_dir) / "metrics" / "evaluation"
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = []
    sample_index = 0
    for batch in loader:
        if len(batch) == 4:
            atm, ocn, atm_mask, ocn_mask = batch
        else:
            atm, ocn = batch
            atm_mask = torch.ones_like(atm, dtype=torch.bool)
            ocn_mask = torch.ones_like(ocn, dtype=torch.bool)
        if max_samples is not None and sample_index >= max_samples:
            break
        if max_samples is not None and sample_index + atm.shape[0] > max_samples:
            keep = max_samples - sample_index
            atm = atm[:keep]
            ocn = ocn[:keep]
            atm_mask = atm_mask[:keep]
            ocn_mask = ocn_mask[:keep]
        out = model(atm.to(device), ocn.to(device))
        atm_recon = out["atm_recon"].cpu()
        ocn_recon = out["ocn_recon"].cpu()

        for batch_idx in range(atm.shape[0]):
            for channel, variable in enumerate(dataset.atm_vars):
                original = _masked_numpy(_denorm(dataset, "atm", variable, atm[batch_idx, channel]), atm_mask[batch_idx, channel])
                reconstructed = _masked_numpy(
                    _denorm(dataset, "atm", variable, atm_recon[batch_idx, channel]),
                    atm_mask[batch_idx, channel],
                )
                metrics.append(
                    compute_variable_metric(
                        original,
                        reconstructed,
                        domain="atm",
                        variable=variable,
                        sample=sample_index + batch_idx,
                    )
                )
            for channel, variable in enumerate(dataset.ocn_vars):
                original = _masked_numpy(_denorm(dataset, "ocn", variable, ocn[batch_idx, channel]), ocn_mask[batch_idx, channel])
                reconstructed = _masked_numpy(
                    _denorm(dataset, "ocn", variable, ocn_recon[batch_idx, channel]),
                    ocn_mask[batch_idx, channel],
                )
                metrics.append(
                    compute_variable_metric(
                        original,
                        reconstructed,
                        domain="ocn",
                        variable=variable,
                        sample=sample_index + batch_idx,
                    )
                )
        sample_index += atm.shape[0]

    summary = summarize_metrics(metrics)
    (output_dir / "per_sample_variable_metrics.csv").write_text(
        metric_header() + "\n" + "\n".join(metric_to_csv(item) for item in metrics),
        encoding="utf-8",
    )
    (output_dir / "summary_metrics.csv").write_text(summary_to_csv(summary), encoding="utf-8")
    (output_dir / "summary_metrics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", default=str(DEFAULT_AE_MODEL_PATH))
    parser.add_argument("--data", default=str(DEFAULT_PROCESSED_PT))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA_PATH))
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--split", choices=["train", "val", "all"], default="val")
    parser.add_argument("--val-fraction", type=float, default=0.1)
    args = parser.parse_args()
    evaluate(
        args.checkpoint,
        args.data,
        args.metadata,
        args.output_dir,
        args.batch_size,
        args.max_samples,
        split=args.split,
        val_fraction=args.val_fraction,
    )


if __name__ == "__main__":
    main()
