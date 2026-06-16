from __future__ import annotations

from typing import Any

import torch
import torch.nn as nn


class EMAModel:
    def __init__(self, model: nn.Module, decay: float = 0.9999, warmup_steps: int = 1000) -> None:
        self.decay = float(decay)
        self.warmup_steps = int(warmup_steps)
        self.step_count = 0
        self.shadow_params = [param.detach().clone() for param in model.parameters()]
        self._backup: list[torch.Tensor] | None = None

    def _get_decay(self) -> float:
        if self.warmup_steps <= 0:
            return self.decay
        return self.decay * min(float(self.step_count) / float(self.warmup_steps), 1.0)

    @torch.no_grad()
    def update(self, model: nn.Module) -> None:
        decay = self._get_decay()
        self.step_count += 1
        for shadow, param in zip(self.shadow_params, model.parameters()):
            shadow.mul_(decay).add_(param.detach(), alpha=1.0 - decay)

    def state_dict(self) -> dict[str, Any]:
        return {
            "decay": self.decay,
            "warmup_steps": self.warmup_steps,
            "step_count": self.step_count,
            "shadow_params": [param.detach().clone() for param in self.shadow_params],
        }

    def load_state_dict(self, state_dict: dict[str, Any], model: nn.Module) -> None:
        shadow_params = state_dict.get("shadow_params")
        if not isinstance(shadow_params, list):
            raise ValueError("EMA state is missing shadow_params.")
        model_params = list(model.parameters())
        if len(shadow_params) != len(model_params):
            raise ValueError(f"EMA parameter count mismatch: checkpoint={len(shadow_params)} model={len(model_params)}")
        self.decay = float(state_dict.get("decay", self.decay))
        self.warmup_steps = int(state_dict.get("warmup_steps", self.warmup_steps))
        self.step_count = int(state_dict.get("step_count", self.step_count))
        self.shadow_params = [
            shadow.detach().to(device=param.device, dtype=param.dtype).clone()
            for shadow, param in zip(shadow_params, model_params)
        ]
        self._backup = None

    def reset_to_model(self, model: nn.Module) -> None:
        self.shadow_params = [param.detach().clone() for param in model.parameters()]
        self._backup = None

    @torch.no_grad()
    def apply_shadow(self, model: nn.Module) -> None:
        self._backup = [param.detach().clone() for param in model.parameters()]
        for shadow, param in zip(self.shadow_params, model.parameters()):
            param.data.copy_(shadow)

    @torch.no_grad()
    def restore(self, model: nn.Module) -> None:
        if self._backup is None:
            return
        for backup, param in zip(self._backup, model.parameters()):
            param.data.copy_(backup)
        self._backup = None
