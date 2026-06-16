from __future__ import annotations

import math
from typing import Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F

from .attention import CrossAttention2d
from .blocks import ConvNeXtBlock, Downsample2d, LayerNorm2d, ResBlock2d, Upsample2d


ATM_SIZE = (90, 144)
OCN_SIZE = (200, 360)
COUPLING_GRID = (50, 90)
LATENT_GRID = (25, 45)


class CrossFusionBlock(nn.Module):
    def __init__(self, channels: int, heads: int = 4, attention_backend: str = "auto") -> None:
        super().__init__()
        self.atm_from_ocn = CrossAttention2d(channels, heads=heads, backend=attention_backend)
        self.ocn_from_atm = CrossAttention2d(channels, heads=heads, backend=attention_backend)
        self.atm_ffn = ConvNeXtBlock(channels)
        self.ocn_ffn = ConvNeXtBlock(channels)

    def forward(self, atm: torch.Tensor, ocn: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # Bidirectional cross-attention is the explicit atmosphere-ocean coupling step.
        atm = self.atm_from_ocn(atm, ocn)
        ocn = self.ocn_from_atm(ocn, atm)
        return self.atm_ffn(atm), self.ocn_ffn(ocn)


class CNNFusionBlock(nn.Module):
    def __init__(self, channels: int, depth: int = 2) -> None:
        super().__init__()
        fused_channels = channels * 2
        self.net = nn.Sequential(
            LayerNorm2d(fused_channels),
            nn.Conv2d(fused_channels, fused_channels, kernel_size=3, padding=1),
            nn.SiLU(),
            *[ResBlock2d(fused_channels) for _ in range(max(1, int(depth)))],
            LayerNorm2d(fused_channels),
            nn.Conv2d(fused_channels, fused_channels, kernel_size=3, padding=1),
        )

    def forward(self, atm: torch.Tensor, ocn: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        delta_atm, delta_ocn = self.net(torch.cat((atm, ocn), dim=1)).chunk(2, dim=1)
        return atm + delta_atm, ocn + delta_ocn


class DomainStem(nn.Module):
    """A shallow per-domain stem before both domains are placed on the latent grid."""

    def __init__(
        self,
        in_channels: int,
        width: int,
        down_kernel: int,
        extra_blocks: int = 0,
        res_blocks_per_level: int = 2,
    ) -> None:
        super().__init__()
        depth = max(1, int(res_blocks_per_level))
        blocks: list[nn.Module] = [
            nn.Conv2d(in_channels, width // 2, kernel_size=3, padding=1),
            *[ResBlock2d(width // 2) for _ in range(depth)],
            Downsample2d(width // 2, width, kernel_size=down_kernel),
            *[ResBlock2d(width) for _ in range(depth)],
        ]
        blocks.extend(ResBlock2d(width) for _ in range(extra_blocks))
        self.net = nn.Sequential(*blocks)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class CoupledEncoder(nn.Module):
    def __init__(
        self,
        in_atm_channels: int = 4,
        in_ocn_channels: int = 4,
        width: int = 56,
        latent_channels: int = 36,
        cross_depth: int = 4,
        heads: int = 4,
        attention_backend: str = "auto",
        fusion_mode: str = "cross_attention",
        cnn_fusion_depth: int = 2,
    ) -> None:
        super().__init__()
        pre_depth = max(1, cross_depth // 2)
        post_depth = max(1, cross_depth - pre_depth)
        fusion_mode = str(fusion_mode).lower()
        if fusion_mode not in {"cross_attention", "attention", "cnn"}:
            raise ValueError(f"Unsupported fusion_mode={fusion_mode!r}; expected 'cross_attention' or 'cnn'.")
        fusion_block = (
            (lambda: CNNFusionBlock(width, depth=cnn_fusion_depth))
            if fusion_mode == "cnn"
            else (lambda: CrossFusionBlock(width, heads=heads, attention_backend=attention_backend))
        )
        self.fusion_mode = fusion_mode
        self.atm_stem = DomainStem(in_atm_channels, width, down_kernel=2)
        self.ocn_stem = DomainStem(in_ocn_channels, width, down_kernel=3, extra_blocks=1)
        self.ocn_to_coupling = nn.Sequential(
            Downsample2d(width, width, kernel_size=3),
            ResBlock2d(width),
            ResBlock2d(width),
        )
        self.atm_refine = nn.Sequential(ResBlock2d(width), ResBlock2d(width))
        self.ocn_refine = nn.Sequential(ResBlock2d(width), ResBlock2d(width))
        self.pre_down_coupling = nn.ModuleList(
            [fusion_block() for _ in range(pre_depth)]
        )
        self.final_atm_down = nn.Sequential(
            Downsample2d(width, width, kernel_size=2),
            ResBlock2d(width),
            ResBlock2d(width),
        )
        self.final_ocn_down = nn.Sequential(
            Downsample2d(width, width, kernel_size=2),
            ResBlock2d(width),
            ResBlock2d(width),
        )
        self.post_down_coupling = nn.ModuleList(
            [fusion_block() for _ in range(post_depth)]
        )
        self.fuse = nn.Sequential(
            LayerNorm2d(width * 2),
            nn.Conv2d(width * 2, width * 2, 3, padding=1),
            nn.SiLU(),
            ResBlock2d(width * 2),
            ResBlock2d(width * 2),
            nn.Conv2d(width * 2, latent_channels, 3, padding=1),
        )

    def forward(self, atm: torch.Tensor, ocn: torch.Tensor) -> torch.Tensor:
        atm = self.atm_stem(atm)
        ocn = self.ocn_stem(ocn)
        atm = self.atm_refine(F.interpolate(atm, size=COUPLING_GRID, mode="bilinear", align_corners=False))
        ocn = self.ocn_refine(self.ocn_to_coupling(ocn))
        for block in self.pre_down_coupling:
            atm, ocn = block(atm, ocn)
        atm = self.final_atm_down(atm)
        ocn = self.final_ocn_down(ocn)
        for block in self.post_down_coupling:
            atm, ocn = block(atm, ocn)
        return self.fuse(torch.cat((atm, ocn), dim=1))


def _decoder_stage_targets(target_size: tuple[int, int], max_ratio: float = 1.5) -> list[tuple[int, int]]:
    targets: list[tuple[int, int]] = []
    current = LATENT_GRID
    while current != target_size:
        next_size = (
            min(target_size[0], max(current[0] + 1, int(current[0] * max_ratio))),
            min(target_size[1], max(current[1] + 1, int(current[1] * max_ratio))),
        )
        if next_size == current:
            next_size = target_size
        targets.append(next_size)
        current = next_size
    return targets


def _expand_stage_depths(stage_depth: int | Sequence[int], num_stages: int) -> list[int]:
    if isinstance(stage_depth, int):
        return [int(stage_depth)] * num_stages
    depths = [int(item) for item in stage_depth]
    if not depths:
        raise ValueError("decoder_stage_depth must contain at least one value")
    if len(depths) < num_stages:
        depths.extend([depths[-1]] * (num_stages - len(depths)))
    return depths[:num_stages]


class DomainDecoder(nn.Module):
    def __init__(
        self,
        latent_channels: int,
        out_channels: int,
        target_size: tuple[int, int],
        width: int,
        middle_depth: int = 4,
        stage_depth: int | Sequence[int] = 2,
        max_resize_ratio: float = 1.5,
    ) -> None:
        super().__init__()
        if width % 8 != 0:
            raise ValueError("decoder width must be divisible by 8")
        self.target_size = target_size
        self.width = int(width)
        self.middle_depth = int(middle_depth)
        self.stage_targets = _decoder_stage_targets(target_size, max_ratio=max_resize_ratio)
        self.stage_depths = _expand_stage_depths(stage_depth, len(self.stage_targets))
        self.stage_depth = self.stage_depths
        self.in_proj = nn.Conv2d(latent_channels, width, 3, padding=1)
        self.middle = nn.Sequential(*[ResBlock2d(width) for _ in range(self.middle_depth)])
        self.ups = nn.ModuleList()
        self.blocks = nn.ModuleList()
        in_width = width
        for index, depth in enumerate(self.stage_depths):
            out_width = max(width // (2 ** (index + 1)), width // 8)
            self.ups.append(Upsample2d(in_width, out_width))
            self.blocks.append(nn.Sequential(*[ResBlock2d(out_width) for _ in range(depth)]))
            in_width = out_width
        self.head = nn.Sequential(
            LayerNorm2d(in_width),
            nn.SiLU(),
            nn.Conv2d(in_width, out_channels, 3, padding=1),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        x = self.middle(self.in_proj(z))
        for target, up, block in zip(self.stage_targets, self.ups, self.blocks):
            x = block(up(x, size=target))
        return self.head(x)


class CoupledAE(nn.Module):
    """Deterministic coupled atmosphere-ocean autoencoder for CM2 LDA."""

    def __init__(
        self,
        latent_dim: int | None = None,
        latent_channels: int | None = None,
        in_atm_channels: int = 4,
        in_ocn_channels: int = 4,
        width: int = 56,
        decoder_width: int = 160,
        decoder_middle_depth: int = 4,
        decoder_stage_depth: int | Sequence[int] = 2,
        atm_decoder_stage_depth: int | Sequence[int] | None = None,
        ocn_decoder_stage_depth: int | Sequence[int] | None = None,
        cross_depth: int = 4,
        heads: int = 4,
        attention_backend: str = "auto",
        fusion_mode: str = "cross_attention",
        cnn_fusion_depth: int = 2,
    ) -> None:
        super().__init__()
        grid_size = LATENT_GRID[0] * LATENT_GRID[1]
        if latent_channels is None:
            latent_channels = max(1, latent_dim // grid_size) if latent_dim and latent_dim % grid_size == 0 else 36
        self.latent_channels = int(latent_channels)
        self.in_atm_channels = int(in_atm_channels)
        self.in_ocn_channels = int(in_ocn_channels)
        self.width = int(width)
        self.decoder_width = int(decoder_width)
        self.decoder_middle_depth = int(decoder_middle_depth)
        self.decoder_stage_depth = list(decoder_stage_depth) if not isinstance(decoder_stage_depth, int) else int(decoder_stage_depth)
        self.atm_decoder_stage_depth = (
            self.decoder_stage_depth
            if atm_decoder_stage_depth is None
            else list(atm_decoder_stage_depth) if not isinstance(atm_decoder_stage_depth, int) else int(atm_decoder_stage_depth)
        )
        self.ocn_decoder_stage_depth = (
            self.decoder_stage_depth
            if ocn_decoder_stage_depth is None
            else list(ocn_decoder_stage_depth) if not isinstance(ocn_decoder_stage_depth, int) else int(ocn_decoder_stage_depth)
        )
        self.cross_depth = int(cross_depth)
        self.heads = int(heads)
        self.attention_backend = str(attention_backend)
        self.fusion_mode = str(fusion_mode)
        self.cnn_fusion_depth = int(cnn_fusion_depth)
        self.latent_shape = (self.latent_channels, *LATENT_GRID)
        self.latent_dim = math.prod(self.latent_shape)
        self.encoder = CoupledEncoder(
            in_atm_channels=in_atm_channels,
            in_ocn_channels=in_ocn_channels,
            width=width,
            latent_channels=self.latent_channels,
            cross_depth=cross_depth,
            heads=heads,
            attention_backend=attention_backend,
            fusion_mode=fusion_mode,
            cnn_fusion_depth=cnn_fusion_depth,
        )
        self.atm_decoder = DomainDecoder(
            self.latent_channels,
            in_atm_channels,
            ATM_SIZE,
            width=self.decoder_width,
            middle_depth=self.decoder_middle_depth,
            stage_depth=self.atm_decoder_stage_depth,
        )
        self.ocn_decoder = DomainDecoder(
            self.latent_channels,
            in_ocn_channels,
            OCN_SIZE,
            width=self.decoder_width,
            middle_depth=self.decoder_middle_depth,
            stage_depth=self.ocn_decoder_stage_depth,
        )

    def _reshape_latent(self, z: torch.Tensor) -> torch.Tensor:
        if z.dim() == 2:
            return z.view(z.shape[0], *self.latent_shape)
        if z.dim() == 3 and z.shape[-1] == 1:
            return z.squeeze(-1).view(z.shape[0], *self.latent_shape)
        return z

    def encode(self, atm: torch.Tensor, ocn: torch.Tensor) -> torch.Tensor:
        return self.encoder(atm, ocn)

    def decode(self, z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        z = self._reshape_latent(z)
        return self.atm_decoder(z), self.ocn_decoder(z)

    def forward(
        self,
        atm: torch.Tensor,
        ocn: torch.Tensor,
        corruption_sigma: float = 0.0,
        single_pass_denoising: bool = False,
    ) -> dict[str, torch.Tensor]:
        latent = self.encode(atm, ocn)
        if single_pass_denoising and corruption_sigma > 0.0:
            decode_latent = latent + float(corruption_sigma) * torch.randn_like(latent)
        else:
            decode_latent = latent
        atm_recon, ocn_recon = self.decode(decode_latent)
        out = {
            "atm_recon": atm_recon,
            "ocn_recon": ocn_recon,
            "latent": latent,
        }
        if corruption_sigma > 0.0 and not single_pass_denoising:
            noisy_latent = latent + float(corruption_sigma) * torch.randn_like(latent)
            atm_corrupt, ocn_corrupt = self.decode(noisy_latent)
            out["atm_corrupt_recon"] = atm_corrupt
            out["ocn_corrupt_recon"] = ocn_corrupt
        return out


def model_config_from_instance(model: CoupledAE) -> dict:
    atm_elements = model.in_atm_channels * math.prod(ATM_SIZE)
    ocn_elements = model.in_ocn_channels * math.prod(OCN_SIZE)
    input_elements = atm_elements + ocn_elements
    return {
        "model_type": "coupled_ae",
        "latent_dim": model.latent_dim,
        "latent_channels": model.latent_channels,
        "in_atm_channels": model.in_atm_channels,
        "in_ocn_channels": model.in_ocn_channels,
        "width": model.width,
        "decoder_width": model.decoder_width,
        "decoder_middle_depth": model.decoder_middle_depth,
        "decoder_stage_depth": model.decoder_stage_depth,
        "atm_decoder_stage_depth": model.atm_decoder_stage_depth,
        "ocn_decoder_stage_depth": model.ocn_decoder_stage_depth,
        "cross_depth": model.cross_depth,
        "heads": model.heads,
        "attention_backend": model.attention_backend,
        "fusion_mode": model.fusion_mode,
        "cnn_fusion_depth": model.cnn_fusion_depth,
        "latent_shape": list(model.latent_shape),
        "atm_size": list(ATM_SIZE),
        "ocn_size": list(OCN_SIZE),
        "latent_grid": list(LATENT_GRID),
        "coupling_grid": list(COUPLING_GRID),
        "compression": {
            "input_elements": input_elements,
            "latent_elements": model.latent_dim,
            "overall_element_ratio": input_elements / model.latent_dim,
            "atm_element_ratio_to_latent": atm_elements / model.latent_dim,
            "ocn_element_ratio_to_latent": ocn_elements / model.latent_dim,
            "atm_spatial_area_ratio_to_latent_grid": math.prod(ATM_SIZE) / math.prod(LATENT_GRID),
            "ocn_spatial_area_ratio_to_latent_grid": math.prod(OCN_SIZE) / math.prod(LATENT_GRID),
            "stem_downsample_steps": 2,
            "domain_downsample_steps_before_latent": 2,
            "max_decoder_resize_ratio_per_stage": 1.5,
        },
    }
