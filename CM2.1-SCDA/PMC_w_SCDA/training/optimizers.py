from __future__ import annotations

from typing import Any

import torch
import torch.nn as nn
from torch.optim import Optimizer
from torch.optim.lr_scheduler import (
    CosineAnnealingLR,
    CosineAnnealingWarmRestarts,
    OneCycleLR,
    ReduceLROnPlateau,
    StepLR,
)


class LAMB(Optimizer):
    def __init__(
        self,
        params,
        lr: float = 1e-3,
        betas: tuple[float, float] = (0.9, 0.999),
        eps: float = 1e-6,
        weight_decay: float = 0.01,
    ) -> None:
        super().__init__(params, dict(lr=lr, betas=betas, eps=eps, weight_decay=weight_decay))

    @torch.no_grad()
    def step(self, closure=None):
        loss = None
        if closure is not None:
            with torch.enable_grad():
                loss = closure()
        for group in self.param_groups:
            for p in group["params"]:
                if p.grad is None:
                    continue
                grad = p.grad
                if grad.is_sparse:
                    raise RuntimeError("LAMB does not support sparse gradients")
                state = self.state[p]
                if len(state) == 0:
                    state["step"] = 0
                    state["exp_avg"] = torch.zeros_like(p)
                    state["exp_avg_sq"] = torch.zeros_like(p)
                exp_avg, exp_avg_sq = state["exp_avg"], state["exp_avg_sq"]
                beta1, beta2 = group["betas"]
                state["step"] += 1
                exp_avg.mul_(beta1).add_(grad, alpha=1.0 - beta1)
                exp_avg_sq.mul_(beta2).addcmul_(grad, grad, value=1.0 - beta2)
                step = (exp_avg / (1.0 - beta1 ** state["step"])) / (
                    (exp_avg_sq / (1.0 - beta2 ** state["step"])).sqrt() + group["eps"]
                )
                if group["weight_decay"] != 0:
                    step.add_(p, alpha=group["weight_decay"])
                trust_ratio = p.norm(2).clamp(min=1e-6) / step.norm(2).clamp(min=1e-6)
                p.add_(step, alpha=-group["lr"] * trust_ratio)
        return loss


class AdEMAMix(Optimizer):
    def __init__(
        self,
        params,
        lr: float = 1e-3,
        betas: tuple[float, float, float] = (0.9, 0.999, 0.9999),
        alpha: float = 5.0,
        eps: float = 1e-8,
        weight_decay: float = 0.0,
    ) -> None:
        super().__init__(params, dict(lr=lr, betas=betas, alpha=alpha, eps=eps, weight_decay=weight_decay))

    @torch.no_grad()
    def step(self, closure=None):
        loss = None
        if closure is not None:
            with torch.enable_grad():
                loss = closure()
        for group in self.param_groups:
            beta1, beta2, beta3 = group["betas"]
            for p in group["params"]:
                if p.grad is None:
                    continue
                grad = p.grad
                state = self.state[p]
                if len(state) == 0:
                    state["step"] = 0
                    state["m_fast"] = torch.zeros_like(p)
                    state["m_slow"] = torch.zeros_like(p)
                    state["v"] = torch.zeros_like(p)
                state["step"] += 1
                if group["weight_decay"] != 0:
                    grad = grad.add(p, alpha=group["weight_decay"])
                state["m_fast"].mul_(beta1).add_(grad, alpha=1.0 - beta1)
                state["m_slow"].mul_(beta3).add_(grad, alpha=1.0 - beta3)
                state["v"].mul_(beta2).addcmul_(grad, grad, value=1.0 - beta2)
                m_hat = (
                    state["m_fast"] / (1.0 - beta1 ** state["step"])
                    + group["alpha"] * state["m_slow"] / (1.0 - beta3 ** state["step"])
                ) / (1.0 + group["alpha"])
                v_hat = state["v"] / (1.0 - beta2 ** state["step"])
                p.addcdiv_(m_hat, v_hat.sqrt().add_(group["eps"]), value=-group["lr"])
        return loss


def build_optimizer(model: nn.Module, cfg: dict[str, Any]) -> Optimizer:
    name = str(cfg.get("name", "adamw")).lower()
    lr = float(cfg.get("lr", 2e-4))
    weight_decay = float(cfg.get("weight_decay", 1e-4))
    betas = tuple(float(item) for item in cfg.get("betas", [0.9, 0.999]))
    eps = float(cfg.get("eps", 1e-8))
    params = list(model.parameters())
    if name == "adamw":
        kwargs: dict[str, Any] = {"lr": lr, "weight_decay": weight_decay, "betas": betas, "eps": eps}
        if "foreach" in cfg:
            kwargs["foreach"] = bool(cfg["foreach"])
        if "capturable" in cfg:
            kwargs["capturable"] = bool(cfg["capturable"])
        if bool(cfg.get("fused", False)) and any(param.is_cuda for param in params):
            kwargs["fused"] = True
        return torch.optim.AdamW(params, **kwargs)
    if name == "adam":
        kwargs = {"lr": lr, "weight_decay": weight_decay, "betas": betas, "eps": eps}
        if "foreach" in cfg:
            kwargs["foreach"] = bool(cfg["foreach"])
        return torch.optim.Adam(params, **kwargs)
    if name == "sgd":
        return torch.optim.SGD(
            params,
            lr=lr,
            weight_decay=weight_decay,
            momentum=float(cfg.get("momentum", 0.9)),
        )
    if name == "lamb":
        return LAMB(params, lr=lr, weight_decay=weight_decay, betas=betas, eps=eps)
    if name == "ademamix":
        betas3 = tuple(float(item) for item in cfg.get("betas", [0.9, 0.999, 0.9999]))
        return AdEMAMix(
            params,
            lr=lr,
            weight_decay=weight_decay,
            betas=betas3,
            alpha=float(cfg.get("alpha", 5.0)),
            eps=eps,
        )
    raise ValueError(f"Unknown optimizer: {name}")


def build_scheduler(optimizer: Optimizer, cfg: dict[str, Any], total_epochs: int, steps_per_epoch: int):
    name = str(cfg.get("name", "plateau")).lower()
    if name == "none":
        return None
    if name == "cosine":
        return CosineAnnealingLR(
            optimizer,
            T_max=max(int(cfg.get("step_size", total_epochs)), 1),
            eta_min=float(cfg.get("min_lr", cfg.get("eta_min", 1e-7))),
        )
    if name == "cosine_warm_restarts":
        return CosineAnnealingWarmRestarts(
            optimizer,
            T_0=int(cfg.get("T_0", 50)),
            T_mult=int(cfg.get("T_mult", 2)),
        )
    if name == "step":
        return StepLR(
            optimizer,
            step_size=int(cfg.get("step_size", 50)),
            gamma=float(cfg.get("gamma", 0.5)),
        )
    if name == "plateau":
        return ReduceLROnPlateau(
            optimizer,
            mode="min",
            patience=int(cfg.get("patience", 5)),
            factor=float(cfg.get("factor", 0.5)),
        )
    if name == "onecycle":
        return OneCycleLR(
            optimizer,
            max_lr=float(cfg.get("max_lr", optimizer.defaults.get("lr", 2e-4))),
            total_steps=max(total_epochs * max(steps_per_epoch, 1), 1),
        )
    raise ValueError(f"Unknown scheduler: {name}")
