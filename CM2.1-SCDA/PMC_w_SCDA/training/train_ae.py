from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import torch
import torch.distributed as dist

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from training.trainer import CoupledAETrainer
from utils.data import build_loaders
from utils.logger import configure_logger, get_logger
from utils.paths import DEFAULT_METADATA_PATH, DEFAULT_PROCESSED_PT, ensure_project_dirs


def _load_yaml(path: str | Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required when Hydra is unavailable.") from exc
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _to_dict(cfg: Any) -> dict[str, Any]:
    if isinstance(cfg, dict):
        return cfg
    try:
        from omegaconf import OmegaConf

        return OmegaConf.to_container(cfg, resolve=True)
    except Exception:
        return dict(cfg)


def _init_distributed(cfg: dict[str, Any]) -> tuple[str, bool, int, int, int]:
    backend = str(cfg.get("distributed", {}).get("backend", "none")).lower()
    if backend in {"none", "false", ""}:
        return "none", False, 0, 0, 1
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    if world_size <= 1:
        return backend, False, 0, 0, 1
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))
    rank = int(os.environ.get("RANK", str(local_rank)))
    if torch.cuda.is_available():
        torch.cuda.set_device(local_rank)
    if backend == "ddp" and not dist.is_initialized():
        dist.init_process_group(backend="nccl" if torch.cuda.is_available() else "gloo")
        rank = dist.get_rank()
    return backend, True, local_rank, rank, world_size


def run_training(cfg: dict[str, Any]) -> Path:
    ensure_project_dirs()
    backend, distributed, local_rank, rank, world_size = _init_distributed(cfg)
    device = torch.device(f"cuda:{local_rank}" if torch.cuda.is_available() else "cpu")
    log_cfg = cfg.get("logging", {})
    configure_logger(
        run_name=str(log_cfg.get("run_name", "cm2_ae")),
        output_dir=log_cfg.get("output_dir"),
        console_level=log_cfg.get("console_level", "INFO"),
        file_level=log_cfg.get("file_level", "DEBUG"),
    )
    logger = get_logger()
    logger.info(
        "Starting CM2-LDA AE training: backend=%s distributed=%s rank=%d local_rank=%d world_size=%d device=%s",
        backend,
        distributed,
        rank,
        local_rank,
        world_size,
        device,
    )

    data_cfg = cfg.get("data", {})
    logger.info(
        "Building data loaders: data.path=%s metadata=%s batch_size=%s val_fraction=%s num_workers=%s prefetch_factor=%s persistent_workers=%s pin_memory=%s timeout=%s drop_last=%s shuffle_strategy=%s block_size=%s block_within=%s mmap_frames=%s in_order=%s worker_torch_threads=%s max_samples=%s",
        data_cfg.get("path", str(DEFAULT_PROCESSED_PT)),
        data_cfg.get("metadata_path", str(DEFAULT_METADATA_PATH)),
        int(cfg.get("training", {}).get("batch_size", 2)),
        data_cfg.get("val_fraction", 0.1),
        data_cfg.get("num_workers", 0),
        data_cfg.get("prefetch_factor", 2),
        data_cfg.get("persistent_workers"),
        data_cfg.get("pin_memory", True),
        data_cfg.get("timeout", 0),
        data_cfg.get("drop_last", True),
        data_cfg.get("shuffle_strategy", "block"),
        data_cfg.get("block_shuffle_size", 512),
        data_cfg.get("block_shuffle_within", False),
        data_cfg.get("mmap_frames", False),
        data_cfg.get("in_order", True),
        data_cfg.get("worker_torch_threads", 1),
        data_cfg.get("max_samples"),
    )
    batch_size = int(cfg.get("training", {}).get("batch_size", 2))
    num_workers = int(data_cfg.get("num_workers", 0))
    prefetch_factor = int(data_cfg.get("prefetch_factor", 2))
    block_size = int(data_cfg.get("block_shuffle_size", 512))
    if str(data_cfg.get("shuffle_strategy", "")).lower() in {"block", "block_shuffle", "chunk", "chunk_shuffle"}:
        rank_chunk = block_size // max(world_size, 1)
        logger.info(
            "Block shuffle geometry: block_samples=%d rank_chunk_samples=%d rank_chunk_batches=%d prefetch_batches_per_rank=%d",
            block_size,
            rank_chunk,
            rank_chunk // max(batch_size, 1),
            num_workers * prefetch_factor if num_workers > 0 else 0,
        )
    train_loader, val_loader, dataset = build_loaders(
        data_cfg.get("path", str(DEFAULT_PROCESSED_PT)),
        data_cfg.get("metadata_path", str(DEFAULT_METADATA_PATH)),
        batch_size=batch_size,
        val_fraction=float(data_cfg.get("val_fraction", 0.1)),
        num_workers=int(data_cfg.get("num_workers", 0)),
        max_samples=data_cfg.get("max_samples"),
        distributed=distributed,
        rank=rank,
        world_size=world_size,
        seed=int(data_cfg.get("seed", 42)),
        shuffle_split=bool(data_cfg.get("shuffle_split", False)),
        shuffle_train=bool(data_cfg.get("shuffle_train", False)),
        shuffle_strategy=str(data_cfg.get("shuffle_strategy", "block")),
        block_shuffle_size=int(data_cfg.get("block_shuffle_size", 512)),
        block_shuffle_within=bool(data_cfg.get("block_shuffle_within", False)),
        mmap_frames=bool(data_cfg.get("mmap_frames", False)),
        in_order=bool(data_cfg.get("in_order", True)),
        worker_torch_threads=int(data_cfg.get("worker_torch_threads", 1)),
        pin_memory=bool(data_cfg.get("pin_memory", True)),
        prefetch_factor=int(data_cfg.get("prefetch_factor", 2)),
        persistent_workers=data_cfg.get("persistent_workers"),
        timeout=float(data_cfg.get("timeout", 0.0)),
        drop_last=bool(data_cfg.get("drop_last", True)),
    )
    logger.info(
        "Data loaders ready: train_batches=%d val_batches=%d dataset=%s train_samples=%s val_samples=%s train_sampler=%s",
        len(train_loader),
        len(val_loader),
        type(dataset).__name__,
        len(getattr(train_loader, "dataset", [])),
        len(getattr(val_loader, "dataset", [])),
        type(getattr(train_loader, "sampler", None)).__name__,
    )
    trainer = CoupledAETrainer(
        cfg=cfg,
        train_loader=train_loader,
        val_loader=val_loader,
        dataset=dataset,
        device=device,
        backend=backend,
        distributed=distributed,
        local_rank=local_rank,
    )
    resume_path = cfg.get("checkpoint", {}).get("resume_from")
    if resume_path and Path(resume_path).exists():
        logger.info("Resume requested from checkpoint.resume_from=%s", resume_path)
        trainer.resume_from_checkpoint(resume_path)
    elif resume_path:
        logger.warning("Resume checkpoint does not exist, starting fresh: %s", resume_path)
    output = trainer.train()
    if backend == "ddp" and distributed and dist.is_initialized():
        dist.destroy_process_group()
    return output


def _main_without_hydra() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "train_ae.yaml"))
    parser.add_argument("--resume-from", default=None)
    parser.add_argument("--local_rank", type=int, default=None)
    args = parser.parse_args()
    cfg = _load_yaml(args.config)
    if args.resume_from:
        cfg.setdefault("checkpoint", {})["resume_from"] = args.resume_from
    run_training(cfg)


try:
    import hydra
    from omegaconf import DictConfig

    @hydra.main(version_base="1.3", config_path="../configs", config_name="train_ae")
    def main(cfg: DictConfig) -> None:
        run_training(_to_dict(cfg))

except Exception:
    main = _main_without_hydra


if __name__ == "__main__":
    if "--config" in sys.argv:
        _main_without_hydra()
    else:
        main()
