from __future__ import annotations

import json
import math
from functools import partial
from inspect import signature
from pathlib import Path
from typing import Any, Iterator, Sequence

import torch
from torch.utils.data import DataLoader, Dataset, DistributedSampler, Sampler, random_split


RUNTIME_ATM_VARS = ["TEMP", "UCOMP", "VCOMP", "PS"]
RUNTIME_OCN_VARS = ["SST", "U_SURF", "V_SURF", "ETA_T"]
FRAME_MANIFEST = "_frame_manifest.json"
_DATALOADER_PARAMS = set(signature(DataLoader).parameters)


class BlockShuffleSampler(Sampler[int]):
    """Shuffle contiguous frame blocks while preserving mostly sequential IO inside each block."""

    def __init__(
        self,
        dataset: Dataset,
        *,
        num_replicas: int = 1,
        rank: int = 0,
        shuffle: bool = True,
        seed: int = 42,
        block_size: int = 512,
        shuffle_within_block: bool = False,
        drop_last: bool = False,
    ) -> None:
        if num_replicas <= 0:
            raise ValueError("num_replicas must be positive")
        if rank < 0 or rank >= num_replicas:
            raise ValueError("rank must be in [0, num_replicas)")
        self.dataset = dataset
        self.num_replicas = int(num_replicas)
        self.rank = int(rank)
        self.shuffle = bool(shuffle)
        self.seed = int(seed)
        self.block_size = max(int(block_size), 1)
        self.shuffle_within_block = bool(shuffle_within_block)
        self.drop_last = bool(drop_last)
        self.epoch = 0
        dataset_len = len(self.dataset)
        self.num_samples = dataset_len // self.num_replicas if self.drop_last else math.ceil(dataset_len / self.num_replicas)
        self.total_size = self.num_samples * self.num_replicas

    def __iter__(self) -> Iterator[int]:
        n = len(self.dataset)
        generator = torch.Generator().manual_seed(self.seed + self.epoch)
        blocks = [list(range(start, min(start + self.block_size, n))) for start in range(0, n, self.block_size)]
        if self.shuffle:
            block_order = torch.randperm(len(blocks), generator=generator).tolist()
        else:
            block_order = list(range(len(blocks)))
        indices: list[int] = []
        for block_idx in block_order:
            block = blocks[block_idx]
            if self.shuffle and self.shuffle_within_block and len(block) > 1:
                local_order = torch.randperm(len(block), generator=generator).tolist()
                block = [block[i] for i in local_order]
            if self.num_replicas > 1:
                chunk_size = math.ceil(len(block) / self.num_replicas)
                start = self.rank * chunk_size
                end = min(start + chunk_size, len(block))
                indices.extend(block[start:end])
            else:
                indices.extend(block)
        if len(indices) < self.num_samples and indices:
            padding_size = self.num_samples - len(indices)
            indices += (indices * math.ceil(padding_size / len(indices)))[:padding_size]
        indices = indices[: self.num_samples]
        return iter(indices)

    def __len__(self) -> int:
        return self.num_samples

    def set_epoch(self, epoch: int) -> None:
        self.epoch = int(epoch)


def _set_worker_torch_threads(num_threads: int, _worker_id: int) -> None:
    if num_threads > 0:
        torch.set_num_threads(int(num_threads))


def _load_frame_payload(path: Path, *, mmap: bool = False) -> dict[str, Any]:
    if not mmap:
        return torch.load(path, map_location="cpu", weights_only=True)
    try:
        return torch.load(path, map_location="cpu", weights_only=True, mmap=True)
    except TypeError:
        return torch.load(path, map_location="cpu", weights_only=True)
    except RuntimeError as exc:
        if "mmap" not in str(exc).lower():
            raise
        return torch.load(path, map_location="cpu", weights_only=True)


def _load_metadata(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_normalization_stats(metadata: dict, base_dir: Path | None = None) -> dict[str, dict[str, float]]:
    raw = metadata.get("normalization_stats")
    if isinstance(raw, dict):
        return dict(raw)
    stats_path = metadata.get("normalization_stats_path") or metadata.get("normalization_stats_file")
    if stats_path:
        path = Path(stats_path)
        if not path.is_absolute() and base_dir is not None:
            path = base_dir / path
        if path.exists():
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    return {}


def _reorder_channels(tensor: torch.Tensor, src_vars: Sequence[str], dst_vars: Sequence[str]) -> torch.Tensor:
    if not src_vars or len(src_vars) != tensor.shape[1]:
        return tensor
    lookup = {name.upper(): i for i, name in enumerate(src_vars)}
    indices = []
    for name in dst_vars:
        if name.upper() not in lookup:
            return tensor
        indices.append(lookup[name.upper()])
    return tensor[:, indices]


def compute_stats(tensor: torch.Tensor, var_names: Sequence[str], domain: str) -> dict[str, dict[str, float]]:
    stats: dict[str, dict[str, float]] = {}
    for idx, name in enumerate(var_names):
        values = tensor[:, idx].float()
        values = values[torch.isfinite(values)]
        if values.numel() == 0:
            mean, std = 0.0, 1.0
        else:
            mean = float(values.mean())
            std = float(values.std(unbiased=False).clamp_min(1e-8))
        stats[f"{domain}::{name}"] = {"mean": mean, "std": std}
    return stats


def normalize_domain(tensor: torch.Tensor, vars_: Sequence[str], domain: str, stats: dict[str, dict[str, float]]) -> torch.Tensor:
    out = tensor.clone()
    for idx, name in enumerate(vars_):
        stat = stats[f"{domain}::{name}"]
        out[:, idx] = (out[:, idx] - stat["mean"]) / max(stat["std"], 1e-8)
    return torch.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def _fallback_domain_mask(tensor: torch.Tensor, domain: str) -> torch.Tensor:
    finite = torch.isfinite(tensor)
    if domain == "ocn":
        nonzero_grid = tensor.ne(0).any(dim=0, keepdim=True)
        return finite & nonzero_grid.expand_as(tensor)
    return finite


class CM2FrameDataset(Dataset):
    """One-file-per-frame dataset, modelled after PALM-GDA's frame manifest infra."""

    def __init__(
        self,
        data_dir: str | Path,
        split: str = "train",
        val_fraction: float = 0.1,
        max_samples: int | None = None,
        seed: int = 42,
        shuffle_split: bool = False,
        manifest: dict | None = None,
        normalization_stats: dict[str, dict[str, float]] | None = None,
        indices: Sequence[int] | None = None,
        load_mmap: bool = False,
    ) -> None:
        self.data_dir = Path(data_dir)
        manifest_path = self.data_dir / FRAME_MANIFEST
        if manifest is None and not manifest_path.exists():
            raise FileNotFoundError(f"{FRAME_MANIFEST} not found in {self.data_dir}")
        if manifest is None:
            with open(manifest_path, "r", encoding="utf-8") as f:
                manifest = json.load(f)
        self.manifest = manifest
        self.atm_vars = list(self.manifest["atm_vars"])
        self.ocn_vars = list(self.manifest["ocn_vars"])
        self.normalization_stats = normalization_stats or _load_normalization_stats(self.manifest, self.data_dir)
        if not self.normalization_stats:
            raise KeyError(f"Normalization stats sidecar is missing for {manifest_path}")
        frame_files = list(self.manifest["frame_files"])
        total = len(frame_files)
        if indices is not None:
            selected = [int(item) for item in indices]
        else:
            n_val = max(1, int(total * val_fraction)) if total > 1 else 0
            if split == "all":
                selected = list(range(total))
            else:
                if shuffle_split:
                    generator = torch.Generator().manual_seed(seed)
                    order = torch.randperm(total, generator=generator).tolist()
                else:
                    order = list(range(total))
                if split == "val" and n_val > 0:
                    selected = order[-n_val:]
                else:
                    selected = order[:-n_val] if n_val > 0 else order
        if max_samples is not None:
            selected = selected[:max_samples]
        self.indices = list(selected)
        self.frame_files = [frame_files[i] for i in selected]
        self.load_mmap = bool(load_mmap)
        time_values = self.manifest.get("time_values")
        if time_values is None:
            time_values = range(total)
        self.time_values = [time_values[i] for i in selected]

    def __len__(self) -> int:
        return len(self.frame_files)

    def __getitem__(self, index: int) -> tuple[torch.Tensor, torch.Tensor]:
        payload = _load_frame_payload(self.data_dir / self.frame_files[index], mmap=self.load_mmap)
        atm = payload["atm"].float()
        ocn = payload["ocn"].float()
        atm_mask = payload.get("atm_mask")
        ocn_mask = payload.get("ocn_mask")
        atm_mask = atm_mask.bool() if atm_mask is not None else _fallback_domain_mask(atm, "atm")
        ocn_mask = ocn_mask.bool() if ocn_mask is not None else _fallback_domain_mask(ocn, "ocn")
        return atm, ocn, atm_mask, ocn_mask

    def export_metadata(self, path: str | Path, extra: dict | None = None) -> None:
        path = Path(path)
        stats_path = path.with_name("normalization_stats.json")
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "w", encoding="utf-8") as f:
            json.dump(self.normalization_stats, f, indent=2)
        metadata = {
            "atm_vars": self.atm_vars,
            "ocn_vars": self.ocn_vars,
            "normalization_stats_file": stats_path.name,
            "normalization_stats_path": str(stats_path),
            "normalization": "palm_gda_zscore",
            "handle_nan": True,
            "nan_fill_value": 0.0,
        }
        if extra:
            metadata.update(extra)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)


class CM2DatasetView(Dataset):
    """Index view that preserves CM2 dataset metadata used by analysis code."""

    def __init__(self, dataset: Dataset, indices: Sequence[int]) -> None:
        self.dataset = dataset
        self.indices = [int(item) for item in indices]
        self.atm_vars = list(getattr(dataset, "atm_vars", []))
        self.ocn_vars = list(getattr(dataset, "ocn_vars", []))
        self.normalization_stats = getattr(dataset, "normalization_stats", {})
        base_indices = list(getattr(dataset, "indices", range(len(dataset))))
        self.selected_indices = [base_indices[i] if i < len(base_indices) else i for i in self.indices]
        time_values = getattr(dataset, "time_values", None)
        if time_values is not None:
            self.time_values = [time_values[i] for i in self.indices]

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, index: int):
        return self.dataset[self.indices[index]]

    def export_metadata(self, path: str | Path, extra: dict | None = None) -> None:
        if hasattr(self.dataset, "export_metadata"):
            self.dataset.export_metadata(path, extra=extra)
            return
        raise AttributeError("Wrapped dataset does not support export_metadata")


class CM2TensorDataset(Dataset):
    """Compatibility loader for legacy single-file `{atm, ocn}` datasets."""

    def __init__(self, data_path: str | Path, metadata_path: str | Path | None = None, normalize: bool = True) -> None:
        data_path = Path(data_path)
        metadata_path = Path(metadata_path) if metadata_path is not None else data_path.with_name(f"{data_path.stem}_metadata.json")
        payload = torch.load(data_path, map_location="cpu", weights_only=False)
        if not isinstance(payload, dict) or "atm" not in payload or "ocn" not in payload:
            raise RuntimeError(f"{data_path} must contain a dict with 'atm' and 'ocn' tensors.")
        metadata = _load_metadata(metadata_path)
        atm = torch.as_tensor(payload["atm"], dtype=torch.float32)
        ocn = torch.as_tensor(payload["ocn"], dtype=torch.float32)
        atm = _reorder_channels(atm, metadata.get("atm_vars", []), RUNTIME_ATM_VARS)
        ocn = _reorder_channels(ocn, metadata.get("ocn_vars", []), RUNTIME_OCN_VARS)
        self.atm_vars = RUNTIME_ATM_VARS[: atm.shape[1]]
        self.ocn_vars = RUNTIME_OCN_VARS[: ocn.shape[1]]
        self.normalization_stats = _load_normalization_stats(metadata, metadata_path.parent) or {
            **compute_stats(atm, self.atm_vars, "atm"),
            **compute_stats(ocn, self.ocn_vars, "ocn"),
        }
        if normalize:
            atm = normalize_domain(atm, self.atm_vars, "atm", self.normalization_stats)
            ocn = normalize_domain(ocn, self.ocn_vars, "ocn", self.normalization_stats)
        self.atm_mask = _fallback_domain_mask(atm, "atm")
        self.ocn_mask = _fallback_domain_mask(ocn, "ocn")
        self.atm = torch.nan_to_num(atm, nan=0.0, posinf=0.0, neginf=0.0)
        self.ocn = torch.nan_to_num(ocn, nan=0.0, posinf=0.0, neginf=0.0)
        self.length = min(self.atm.shape[0], self.ocn.shape[0])

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: int) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        return self.atm[index], self.ocn[index], self.atm_mask[index], self.ocn_mask[index]

    def export_metadata(self, path: str | Path, extra: dict | None = None) -> None:
        path = Path(path)
        stats_path = path.with_name("normalization_stats.json")
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "w", encoding="utf-8") as f:
            json.dump(self.normalization_stats, f, indent=2)
        metadata = {
            "atm_vars": self.atm_vars,
            "ocn_vars": self.ocn_vars,
            "normalization_stats_file": stats_path.name,
            "normalization_stats_path": str(stats_path),
            "normalization": "palm_gda_zscore",
            "handle_nan": True,
            "nan_fill_value": 0.0,
        }
        if extra:
            metadata.update(extra)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)


def build_dataset(data_path: str | Path, metadata_path: str | Path | None, split: str, val_fraction: float, max_samples: int | None = None) -> Dataset:
    path = Path(data_path)
    if path.is_dir() and (path / FRAME_MANIFEST).exists():
        return CM2FrameDataset(path, split=split, val_fraction=val_fraction, max_samples=max_samples)
    dataset = CM2TensorDataset(path, metadata_path=metadata_path, normalize=True)
    total = len(dataset)
    if split != "all":
        n_val = max(1, int(total * val_fraction)) if total > 1 else 0
        order = list(range(total))
        indices = order[-n_val:] if split == "val" and n_val > 0 else order[:-n_val] if n_val > 0 else order
        dataset = CM2DatasetView(dataset, indices)
    if max_samples is not None:
        dataset = CM2DatasetView(dataset, range(min(int(max_samples), len(dataset))))
    return dataset


def build_loaders(
    data_path: str | Path,
    metadata_path: str | Path | None,
    batch_size: int,
    val_fraction: float = 0.1,
    num_workers: int = 0,
    max_samples: int | None = None,
    distributed: bool = False,
    rank: int = 0,
    world_size: int = 1,
    seed: int = 42,
    shuffle_split: bool = False,
    shuffle_train: bool = True,
    shuffle_strategy: str = "block",
    block_shuffle_size: int = 512,
    block_shuffle_within: bool = False,
    mmap_frames: bool = False,
    in_order: bool = True,
    worker_torch_threads: int = 1,
    pin_memory: bool = True,
    prefetch_factor: int = 2,
    persistent_workers: bool | None = None,
    timeout: float = 0.0,
    drop_last: bool = True,
) -> tuple[DataLoader, DataLoader, Any]:
    path = Path(data_path)
    if path.is_dir() and (path / FRAME_MANIFEST).exists():
        with open(path / FRAME_MANIFEST, "r", encoding="utf-8") as f:
            manifest = json.load(f)
        normalization_stats = _load_normalization_stats(manifest, path)
        total = len(manifest["frame_files"])
        n_val = max(1, int(total * val_fraction)) if total > 1 else 0
        if shuffle_split:
            generator = torch.Generator().manual_seed(seed)
            order = torch.randperm(total, generator=generator).tolist()
        else:
            order = list(range(total))
        train_indices = order[:-n_val] if n_val > 0 else order
        val_indices = order[-n_val:] if n_val > 0 else []
        if max_samples is not None:
            train_indices = train_indices[:max_samples]
            val_indices = val_indices[:max_samples]
        train_set = CM2FrameDataset(
            path,
            split="train",
            val_fraction=val_fraction,
            seed=seed,
            shuffle_split=shuffle_split,
            manifest=manifest,
            normalization_stats=normalization_stats,
            indices=train_indices,
            load_mmap=bool(mmap_frames),
        )
        val_set = CM2FrameDataset(
            path,
            split="val",
            val_fraction=val_fraction,
            seed=seed,
            shuffle_split=shuffle_split,
            manifest=manifest,
            normalization_stats=normalization_stats,
            indices=val_indices,
            load_mmap=bool(mmap_frames),
        )
        strategy = str(shuffle_strategy or "random").lower()
        if strategy in {"block", "block_shuffle", "chunk", "chunk_shuffle"}:
            train_sampler = BlockShuffleSampler(
                train_set,
                num_replicas=world_size if distributed else 1,
                rank=rank if distributed else 0,
                shuffle=shuffle_train,
                seed=seed,
                block_size=int(block_shuffle_size),
                shuffle_within_block=bool(block_shuffle_within),
                drop_last=bool(drop_last),
            )
        elif distributed:
            train_sampler = DistributedSampler(train_set, num_replicas=world_size, rank=rank, shuffle=shuffle_train, seed=seed)
        else:
            train_sampler = None
        val_sampler = DistributedSampler(val_set, num_replicas=world_size, rank=rank, shuffle=False) if distributed else None
        train_workers = max(0, int(num_workers))
        val_workers = max(0, train_workers // 2)
        persistent = train_workers > 0 if persistent_workers is None else bool(persistent_workers)
        common_train_kwargs = {
            "num_workers": train_workers,
            "pin_memory": pin_memory,
            "drop_last": bool(drop_last),
            "persistent_workers": persistent and train_workers > 0,
            "timeout": float(timeout) if train_workers > 0 else 0,
        }
        common_val_kwargs = {
            "num_workers": val_workers,
            "pin_memory": pin_memory,
            "persistent_workers": persistent and val_workers > 0,
            "timeout": float(timeout) if val_workers > 0 else 0,
        }
        if train_workers > 0:
            common_train_kwargs["prefetch_factor"] = prefetch_factor
            common_train_kwargs["worker_init_fn"] = partial(_set_worker_torch_threads, int(worker_torch_threads))
            if "in_order" in _DATALOADER_PARAMS:
                common_train_kwargs["in_order"] = bool(in_order)
        if val_workers > 0:
            common_val_kwargs["prefetch_factor"] = prefetch_factor
            common_val_kwargs["worker_init_fn"] = partial(_set_worker_torch_threads, int(worker_torch_threads))
            if "in_order" in _DATALOADER_PARAMS:
                common_val_kwargs["in_order"] = bool(in_order)
        train_loader = DataLoader(
            train_set,
            batch_size=batch_size,
            shuffle=train_sampler is None and shuffle_train,
            sampler=train_sampler,
            **common_train_kwargs,
        )
        val_loader = DataLoader(
            val_set,
            batch_size=batch_size,
            shuffle=False,
            sampler=val_sampler,
            **common_val_kwargs,
        )
        return train_loader, val_loader, train_set

    dataset = CM2TensorDataset(path, metadata_path=metadata_path, normalize=True)
    val_len = max(1, int(len(dataset) * val_fraction)) if len(dataset) > 1 else 0
    train_len = len(dataset) - val_len
    if val_len > 0:
        train_set, val_set = random_split(dataset, [train_len, val_len], generator=torch.Generator().manual_seed(42))
    else:
        train_set = dataset
        val_set = dataset
    train_sampler = DistributedSampler(train_set, num_replicas=world_size, rank=rank, shuffle=True) if distributed else None
    val_sampler = DistributedSampler(val_set, num_replicas=world_size, rank=rank, shuffle=False) if distributed else None
    train_loader = DataLoader(
        train_set,
        batch_size=batch_size,
        shuffle=train_sampler is None,
        sampler=train_sampler,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=bool(drop_last),
        timeout=float(timeout) if num_workers > 0 else 0,
    )
    val_loader = DataLoader(
        val_set,
        batch_size=batch_size,
        shuffle=False,
        sampler=val_sampler,
        num_workers=num_workers,
        pin_memory=True,
        timeout=float(timeout) if num_workers > 0 else 0,
    )
    return train_loader, val_loader, dataset
