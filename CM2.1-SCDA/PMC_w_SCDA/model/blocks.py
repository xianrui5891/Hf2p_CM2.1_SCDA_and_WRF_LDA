from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class LayerNorm2d(nn.Module):
    def __init__(self, channels: int, eps: float = 1e-6) -> None:
        super().__init__()
        self.weight = nn.Parameter(torch.ones(channels))
        self.bias = nn.Parameter(torch.zeros(channels))
        self.eps = eps

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        mean = x.mean(dim=1, keepdim=True)
        var = (x - mean).pow(2).mean(dim=1, keepdim=True)
        x = (x - mean) * torch.rsqrt(var + self.eps)
        return x * self.weight[:, None, None] + self.bias[:, None, None]


class ConvNeXtBlock(nn.Module):
    def __init__(self, dim: int, mlp_ratio: int = 4, layer_scale: float = 1e-6) -> None:
        super().__init__()
        self.dwconv = nn.Conv2d(dim, dim, kernel_size=7, padding=3, groups=dim)
        self.norm = nn.LayerNorm(dim, eps=1e-6)
        self.pw1 = nn.Linear(dim, mlp_ratio * dim)
        self.pw2 = nn.Linear(mlp_ratio * dim, dim)
        self.gamma = nn.Parameter(layer_scale * torch.ones(dim)) if layer_scale > 0 else None

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        x = self.dwconv(x)
        x = x.permute(0, 2, 3, 1)
        x = self.norm(x)
        x = self.pw2(F.gelu(self.pw1(x)))
        if self.gamma is not None:
            x = self.gamma * x
        x = x.permute(0, 3, 1, 2).contiguous()
        return residual + x


class ResBlock2d(nn.Module):
    def __init__(self, channels: int, hidden: int | None = None) -> None:
        super().__init__()
        hidden = hidden or channels
        self.net = nn.Sequential(
            LayerNorm2d(channels),
            nn.SiLU(),
            nn.Conv2d(channels, hidden, 3, padding=1),
            LayerNorm2d(hidden),
            nn.SiLU(),
            nn.Conv2d(hidden, channels, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class Downsample2d(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, kernel_size: int = 3) -> None:
        super().__init__()
        padding = 0 if kernel_size % 2 == 0 else kernel_size // 2
        self.op = nn.Sequential(
            LayerNorm2d(in_channels),
            nn.Conv2d(in_channels, out_channels, kernel_size=kernel_size, stride=2, padding=padding),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.op(x)


class Upsample2d(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.proj = nn.Conv2d(in_channels, out_channels, kernel_size=3, padding=1)

    def forward(self, x: torch.Tensor, size: tuple[int, int] | None = None) -> torch.Tensor:
        if size is None:
            x = F.interpolate(x, scale_factor=2, mode="bilinear", align_corners=False)
        else:
            x = F.interpolate(x, size=size, mode="bilinear", align_corners=False)
        return self.proj(x)


# Compatibility names used by the old temp model.
LayerNorm = LayerNorm2d
Block = ConvNeXtBlock
