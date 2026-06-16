from __future__ import annotations

from pathlib import Path
from typing import Any

import torch
from torch.nn.parallel import DistributedDataParallel

from model import model_config_from_instance


def unwrap_model(model: torch.nn.Module) -> torch.nn.Module:
    if isinstance(model, DistributedDataParallel):
        return model.module
    if hasattr(model, "module"):
        return model.module
    if hasattr(model, "_orig_mod"):
        return model._orig_mod
    return model


def save_checkpoint(
    path: str | Path,
    model: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    best_val_loss: float,
    extra: dict[str, Any] | None = None,
) -> None:
    raw_model = unwrap_model(model)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "epoch": int(epoch),
        "best_val_loss": float(best_val_loss),
        "model_state_dict": raw_model.state_dict(),
        "latent_dim": raw_model.latent_dim,
        "model_config": model_config_from_instance(raw_model),
    }
    try:
        payload["optimizer_state_dict"] = optimizer.state_dict()
    except Exception as exc:
        payload["optimizer_state_dict_error"] = f"{type(exc).__name__}: {exc}"
    if extra:
        payload.update(extra)
    torch.save(payload, path)
