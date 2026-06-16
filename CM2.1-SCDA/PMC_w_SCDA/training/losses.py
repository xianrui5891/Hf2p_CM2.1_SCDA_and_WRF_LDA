from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import torch
import torch.nn as nn
import torch.nn.functional as F


_SPECTRAL_GEOMETRY_CACHE: dict[tuple[str, int, int, tuple[tuple[float, float], ...]], dict[str, object]] = {}
_SPATIAL_SUBBAND_MASK_CACHE: dict[tuple[str, int, int, tuple[tuple[float, float], ...], float], torch.Tensor] = {}


def _finite_pair(pred: torch.Tensor, target: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    mask = torch.isfinite(pred) & torch.isfinite(target)
    pred = torch.where(mask, pred, torch.zeros_like(pred))
    target = torch.where(mask, target, torch.zeros_like(target))
    return pred, target


def _valid_mask(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    valid = torch.isfinite(pred) & torch.isfinite(target)
    if mask is not None:
        valid = valid & mask.to(device=pred.device, dtype=torch.bool)
    return valid


def _masked_l1_loss(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    valid = _valid_mask(pred, target, mask)
    diff = torch.where(valid, (pred - target).abs(), torch.zeros_like(pred))
    return diff.sum() / valid.sum().clamp_min(1).to(diff.dtype)


def _masked_mse_loss(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    valid = _valid_mask(pred, target, mask)
    diff = torch.where(valid, (pred - target).square(), torch.zeros_like(pred))
    return diff.sum() / valid.sum().clamp_min(1).to(diff.dtype)


def _mask_pair(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> tuple[torch.Tensor, torch.Tensor]:
    valid = _valid_mask(pred, target, mask)
    return torch.where(valid, pred, torch.zeros_like(pred)), torch.where(valid, target, torch.zeros_like(target))


def gradient_loss_2d(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    pred, target = _finite_pair(pred, target)
    valid = _valid_mask(pred, target, mask)
    dx_p = pred[..., :, 1:] - pred[..., :, :-1]
    dx_t = target[..., :, 1:] - target[..., :, :-1]
    dy_p = pred[..., 1:, :] - pred[..., :-1, :]
    dy_t = target[..., 1:, :] - target[..., :-1, :]
    mask_x = valid[..., :, 1:] & valid[..., :, :-1]
    mask_y = valid[..., 1:, :] & valid[..., :-1, :]
    return _masked_l1_loss(dx_p, dx_t, mask_x) + _masked_l1_loss(dy_p, dy_t, mask_y)


def energy_loss(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    valid = _valid_mask(pred, target, mask)
    pred, target = _mask_pair(pred, target, valid)
    denom = valid.sum(dim=(-2, -1)).clamp_min(1).to(pred.dtype)
    pred_energy = pred.square().sum(dim=(-2, -1)) / denom
    target_energy = target.square().sum(dim=(-2, -1)) / denom
    return F.l1_loss(pred_energy, target_energy)


def spectral_loss_2d(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    pred, target = _mask_pair(pred.float(), target.float(), mask)
    pred_fft = torch.fft.rfft2(pred, norm="ortho").abs()
    target_fft = torch.fft.rfft2(target, norm="ortho").abs()
    return F.l1_loss(torch.log1p(pred_fft), torch.log1p(target_fft))


def _radial_bin_index(height: int, width: int, device: torch.device) -> tuple[torch.Tensor, int]:
    ky = torch.fft.fftfreq(height, device=device)[:, None]
    kx = torch.fft.rfftfreq(width, device=device)[None, :]
    radius = torch.sqrt(ky.square() + kx.square())
    max_radius = radius.max().clamp_min(1e-8)
    bins = torch.clamp((radius / max_radius * (min(height, width) // 2)).long(), min=0)
    return bins.reshape(-1), int(bins.max().item()) + 1


def radial_power_spectrum(x: torch.Tensor) -> torch.Tensor:
    """Return azimuthally averaged 2-D power spectra for `[B, C, H, W]`."""
    x = torch.nan_to_num(x.float(), nan=0.0, posinf=0.0, neginf=0.0)
    _, _, height, width = x.shape
    power = torch.fft.rfft2(x, norm="ortho").abs().square()
    bins, num_bins = _radial_bin_index(height, width, x.device)
    flat = power.reshape(*power.shape[:2], -1)
    spectrum = x.new_zeros((*power.shape[:2], num_bins))
    counts = x.new_zeros(num_bins)
    spectrum.scatter_add_(-1, bins.view(1, 1, -1).expand_as(flat), flat)
    counts.scatter_add_(0, bins, torch.ones_like(bins, dtype=x.dtype))
    return spectrum / counts.clamp_min(1.0)


def radial_spectrum_loss_2d(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
    pred, target = _mask_pair(pred, target, mask)
    pred_spec = radial_power_spectrum(pred)
    target_spec = radial_power_spectrum(target)
    return F.l1_loss(torch.log1p(pred_spec), torch.log1p(target_spec))


def subband_spectrum_loss_2d(
    pred: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor | None = None,
    bands: tuple[tuple[float, float], ...] = ((0.0, 0.15), (0.15, 0.35), (0.35, 1.01)),
    band_weights: tuple[float, ...] = (0.5, 1.0, 1.5),
) -> torch.Tensor:
    """PALM-GDA-inspired spectral band loss over low/mid/high wavenumbers."""
    pred, target = _mask_pair(pred.float(), target.float(), mask)
    _, _, height, width = pred.shape
    ky = torch.fft.fftfreq(height, device=pred.device)[:, None]
    kx = torch.fft.rfftfreq(width, device=pred.device)[None, :]
    radius = torch.sqrt(ky.square() + kx.square())
    radius = radius / radius.max().clamp_min(1e-8)
    pred_power = torch.log1p(torch.fft.rfft2(pred, norm="ortho").abs().square())
    target_power = torch.log1p(torch.fft.rfft2(target, norm="ortho").abs().square())
    total = pred.new_tensor(0.0)
    weight_sum = 0.0
    for (lo, hi), weight in zip(bands, band_weights):
        mask = (radius >= lo) & (radius < hi)
        if mask.any():
            total = total + float(weight) * F.l1_loss(pred_power[..., mask], target_power[..., mask])
            weight_sum += float(weight)
    return total / max(weight_sum, 1e-8)


def _soft_lowpass(radius: torch.Tensor, cutoff: float, half_width: float) -> torch.Tensor:
    if half_width <= 0.0:
        return (radius <= cutoff).to(radius.dtype)
    start = max(0.0, cutoff - half_width)
    end = min(1.0, cutoff + half_width)
    if end <= start:
        return (radius <= cutoff).to(radius.dtype)
    mask = torch.ones_like(radius)
    mask = torch.where(radius >= end, torch.zeros_like(mask), mask)
    transition = (radius > start) & (radius < end)
    if transition.any():
        t = (radius[transition] - start) / (end - start)
        mask[transition] = 0.5 * (1.0 + torch.cos(torch.pi * t))
    return mask


def _get_spatial_subband_masks(
    height: int,
    width: int,
    device: torch.device,
    bands: tuple[tuple[float, float], ...],
    transition_ratio: float,
) -> torch.Tensor:
    key = (str(device), int(height), int(width), bands, float(transition_ratio))
    cached = _SPATIAL_SUBBAND_MASK_CACHE.get(key)
    if cached is not None:
        return cached

    ky = torch.fft.fftfreq(height, device=device)[:, None]
    kx = torch.fft.rfftfreq(width, device=device)[None, :]
    radius = torch.sqrt(ky.square() + kx.square())
    radius = radius / radius.max().clamp_min(1e-8)
    num_bands = len(bands)
    if num_bands <= 1:
        masks = torch.ones(1, height, width // 2 + 1, device=device, dtype=torch.float32)
        _SPATIAL_SUBBAND_MASK_CACHE[key] = masks
        return masks

    edges = [float(bands[0][0]), *[float(hi) for _, hi in bands]]
    internal_lowpasses = []
    for idx in range(1, len(edges) - 1):
        left_width = edges[idx] - edges[idx - 1]
        right_width = edges[idx + 1] - edges[idx]
        half_width = min(float(transition_ratio) * min(left_width, right_width), 0.49 * min(left_width, right_width))
        internal_lowpasses.append(_soft_lowpass(radius, edges[idx], half_width))

    out = []
    for band_idx in range(num_bands):
        if band_idx == 0:
            mask = internal_lowpasses[0]
        elif band_idx == num_bands - 1:
            mask = 1.0 - internal_lowpasses[-1]
        else:
            mask = internal_lowpasses[band_idx] - internal_lowpasses[band_idx - 1]
        out.append(mask.clamp(0.0, 1.0))
    masks = torch.stack(out, dim=0).to(torch.float32).contiguous()
    _SPATIAL_SUBBAND_MASK_CACHE[key] = masks
    return masks


def spatial_subband_loss_2d(
    pred: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor | None = None,
    bands: tuple[tuple[float, float], ...] = ((0.0, 0.25), (0.25, 0.50), (0.50, 0.72), (0.72, 1.01)),
    band_weights: tuple[float, ...] = (0.85, 0.95, 1.30, 1.45),
    transition_ratio: float = 0.12,
    norm: str = "charbonnier",
    relative: bool = False,
    eps: float = 1.0e-3,
) -> torch.Tensor:
    pred, target = _mask_pair(pred.float(), target.float(), mask)
    _, _, height, width = pred.shape
    band_masks = _get_spatial_subband_masks(height, width, pred.device, bands, transition_ratio)
    pred_fft = torch.fft.rfft2(pred, norm="ortho")
    target_fft = torch.fft.rfft2(target, norm="ortho")
    masks = band_masks.to(device=pred.device, dtype=pred_fft.real.dtype).view(1, 1, -1, height, width // 2 + 1)
    pred_bands = torch.fft.irfft2(pred_fft.unsqueeze(2) * masks, s=(height, width), norm="ortho")
    target_bands = torch.fft.irfft2(target_fft.unsqueeze(2) * masks, s=(height, width), norm="ortho")

    diff = pred_bands - target_bands
    if norm == "l1":
        band_error = diff.abs()
    elif norm == "l2":
        band_error = diff.square()
    else:
        band_error = torch.sqrt(diff.square() + float(eps) ** 2)
    band_error = band_error.mean(dim=(1, 3, 4))
    if relative:
        band_scale = target_bands.abs().mean(dim=(1, 3, 4)).clamp_min(float(eps))
        band_error = band_error / band_scale

    weights = pred.new_tensor([float(item) for item in band_weights]).view(1, -1)
    return ((band_error * weights).sum(dim=1) / weights.sum().clamp_min(1e-8)).mean()


def _get_spectral_geometry(
    height: int,
    width: int,
    device: torch.device,
    bands: tuple[tuple[float, float], ...],
) -> dict[str, object]:
    key = (str(device), int(height), int(width), bands)
    cached = _SPECTRAL_GEOMETRY_CACHE.get(key)
    if cached is not None:
        return cached

    ky = torch.fft.fftfreq(height, device=device)[:, None]
    kx = torch.fft.rfftfreq(width, device=device)[None, :]
    radius = torch.sqrt(ky.square() + kx.square())
    norm_radius = radius / radius.max().clamp_min(1e-8)

    num_bins = max(min(height, width) // 2 + 1, 1)
    bins = torch.clamp((norm_radius * (num_bins - 1)).long(), min=0, max=num_bins - 1).reshape(-1)
    counts = torch.zeros(num_bins, device=device, dtype=torch.float32)
    counts.scatter_add_(0, bins, torch.ones_like(bins, dtype=torch.float32))
    band_masks = tuple((norm_radius >= lo) & (norm_radius < hi) for lo, hi in bands)

    cached = {
        "norm_radius": norm_radius,
        "bins": bins,
        "counts": counts.clamp_min(1.0),
        "num_bins": num_bins,
        "band_masks": band_masks,
    }
    _SPECTRAL_GEOMETRY_CACHE[key] = cached
    return cached


def spectral_bundle_loss_2d(
    pred: torch.Tensor,
    target: torch.Tensor,
    *,
    need_spectral: bool,
    need_radial: bool,
    need_subband: bool,
    bands: tuple[tuple[float, float], ...],
    band_weights: tuple[float, ...],
    mask: torch.Tensor | None = None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Compute configured spectral losses from one shared FFT pair."""
    zero = pred.new_tensor(0.0)
    if not need_spectral and not need_radial and not need_subband:
        return zero, zero, zero

    pred, target = _mask_pair(pred.float(), target.float(), mask)
    _, _, height, width = pred.shape
    geometry = _get_spectral_geometry(height, width, pred.device, bands)
    pred_fft = torch.fft.rfft2(pred, norm="ortho").abs()
    target_fft = torch.fft.rfft2(target, norm="ortho").abs()

    spectral = F.l1_loss(torch.log1p(pred_fft), torch.log1p(target_fft)) if need_spectral else zero

    radial = zero
    subband = zero
    if need_radial or need_subband:
        pred_power = pred_fft.square()
        target_power = target_fft.square()

        if need_radial:
            bins = geometry["bins"]
            counts = geometry["counts"]
            num_bins = int(geometry["num_bins"])
            pred_flat = pred_power.reshape(*pred_power.shape[:2], -1)
            target_flat = target_power.reshape(*target_power.shape[:2], -1)
            pred_spec = pred.new_zeros((*pred_power.shape[:2], num_bins))
            target_spec = target.new_zeros((*target_power.shape[:2], num_bins))
            pred_spec.scatter_add_(-1, bins.view(1, 1, -1).expand_as(pred_flat), pred_flat)
            target_spec.scatter_add_(-1, bins.view(1, 1, -1).expand_as(target_flat), target_flat)
            pred_spec = pred_spec / counts
            target_spec = target_spec / counts
            radial = F.l1_loss(torch.log1p(pred_spec), torch.log1p(target_spec))

        if need_subband:
            pred_log_power = torch.log1p(pred_power)
            target_log_power = torch.log1p(target_power)
            total = zero
            weight_sum = 0.0
            for mask, weight in zip(geometry["band_masks"], band_weights):
                total = total + float(weight) * F.l1_loss(pred_log_power[..., mask], target_log_power[..., mask])
                weight_sum += float(weight)
            subband = total / max(weight_sum, 1e-8)

    return spectral, radial, subband


_SCHEDULE_KEYS = ("start_epoch", "ramp_epochs", "end_epoch", "cooldown_epochs", "end_factor")


def get_loss_section(loss_cfg: Mapping[str, Any], section_name: str) -> dict[str, Any]:
    raw = loss_cfg.get(section_name, {})
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raw = {"weight": raw}
    section = dict(raw)
    section.setdefault("enabled", True)
    section.setdefault("start_epoch", 0)
    section.setdefault("ramp_epochs", 0)
    section.setdefault("end_epoch", None)
    section.setdefault("cooldown_epochs", 0)
    section.setdefault("end_factor", 0.0)
    return section


def get_schedule_factor(section_cfg: Mapping[str, Any], epoch: int) -> float:
    if not bool(section_cfg.get("enabled", True)):
        return 0.0
    start_epoch = int(section_cfg.get("start_epoch", 0))
    ramp_epochs = int(section_cfg.get("ramp_epochs", 0))
    if epoch < start_epoch:
        return 0.0
    if ramp_epochs <= 0:
        start_factor = 1.0
    else:
        start_factor = min(float(epoch - start_epoch + 1) / float(ramp_epochs), 1.0)

    end_epoch = section_cfg.get("end_epoch")
    if end_epoch is None:
        return start_factor
    end_epoch = int(end_epoch)
    cooldown_epochs = int(section_cfg.get("cooldown_epochs", 0))
    end_factor_target = float(section_cfg.get("end_factor", 0.0))
    if epoch <= end_epoch:
        end_factor = 1.0
    elif cooldown_epochs <= 0:
        end_factor = end_factor_target
    else:
        progress = min(float(epoch - end_epoch) / float(cooldown_epochs), 1.0)
        end_factor = 1.0 + (end_factor_target - 1.0) * progress
    return start_factor * end_factor


def get_effective_weight(
    section_cfg: Mapping[str, Any],
    epoch: int,
    *,
    key: str = "weight",
    default: float = 0.0,
) -> float:
    base_weight = float(section_cfg.get(key, default))
    if base_weight <= 0.0:
        return 0.0
    return base_weight * get_schedule_factor(section_cfg, epoch)


def _component_section(section: Mapping[str, Any], component: str) -> dict[str, Any]:
    resolved = dict(section)
    for key in _SCHEDULE_KEYS:
        override_key = f"{component}_{key}"
        if override_key in section:
            resolved[key] = section[override_key]
    return resolved


def _subband_bands(section: Mapping[str, Any]) -> tuple[tuple[float, float], ...]:
    edges = section.get("edges")
    if edges is None:
        return ((0.0, 0.15), (0.15, 0.35), (0.35, 1.01))
    values = [float(item) for item in edges]
    if len(values) < 2:
        return ((0.0, 1.01),)
    return tuple((values[i], values[i + 1]) for i in range(len(values) - 1))


def _subband_weights(section: Mapping[str, Any]) -> tuple[float, ...]:
    weights = section.get("weights")
    if weights is None:
        return (0.5, 1.0, 1.5)
    return tuple(float(item) for item in weights)


@dataclass
class LossWeights:
    l1: float = 1.0
    l2: float = 0.25
    gradient: float = 0.08
    energy: float = 0.02
    spectral: float = 0.02
    radial_spectral: float = 0.02
    subband: float = 0.02
    corrupt: float = 0.25


class CoupledAELoss(nn.Module):
    """Small PALM-GDA inspired loss stack for CM2 coupled fields."""

    def __init__(self, weights: LossWeights | None = None, loss_cfg: Mapping[str, Any] | None = None) -> None:
        super().__init__()
        if loss_cfg is None:
            self.weights = weights or LossWeights()
            self._scheduled = False
            self._sections: dict[str, dict[str, Any]] = {}
            self._recon_l1_cfg: dict[str, Any] = {}
            self._recon_l2_cfg: dict[str, Any] = {}
        else:
            self.weights = weights or LossWeights()
            self._scheduled = True
            self._sections = {
                name: get_loss_section(loss_cfg, name)
                for name in ("gradient", "energy", "spectral", "radial_spectral", "subband", "latent_corruption")
            }
            if "weight" not in self._sections["latent_corruption"] and "corrupt" in loss_cfg:
                self._sections["latent_corruption"]["weight"] = loss_cfg["corrupt"]
            recon_cfg = get_loss_section(loss_cfg, "reconstruction")
            if "l1_weight" not in recon_cfg and "l1" in loss_cfg:
                recon_cfg["l1_weight"] = loss_cfg["l1"]
            if "l2_weight" not in recon_cfg and "l2" in loss_cfg:
                recon_cfg["l2_weight"] = loss_cfg["l2"]
            self._recon_l1_cfg = _component_section(recon_cfg, "l1")
            self._recon_l2_cfg = _component_section(recon_cfg, "l2")
            self.apply_epoch_schedule(0)
        subband_cfg = self._sections.get("subband", {}) if self._scheduled else {}
        self.subband_bands = _subband_bands(subband_cfg)
        subband_weights = list(_subband_weights(subband_cfg))
        weight_power = float(subband_cfg.get("weight_power", 0.0))
        if weight_power > 0.0:
            weight_cap = float(subband_cfg.get("weight_cap", 1.0))
            centers = [min(max(0.5 * (lo + hi), 1.0e-8), weight_cap) for lo, hi in self.subband_bands]
            empirical = [center**weight_power for center in centers]
            mean_empirical = sum(empirical) / max(len(empirical), 1)
            subband_weights = [w * e / max(mean_empirical, 1.0e-8) for w, e in zip(subband_weights, empirical)]
        self.subband_weights = tuple(subband_weights)
        self.subband_mode = str(subband_cfg.get("mode", "log_power")).lower()
        self.subband_transition_ratio = float(subband_cfg.get("transition_ratio", 0.12))
        self.subband_norm = str(subband_cfg.get("norm", "charbonnier")).lower()
        self.subband_relative = bool(subband_cfg.get("relative", False))
        self.subband_eps = float(subband_cfg.get("eps", 1.0e-3))

    def apply_epoch_schedule(self, epoch: int) -> None:
        if not self._scheduled:
            return
        self.weights.l1 = get_effective_weight(self._recon_l1_cfg, epoch, key="l1_weight", default=1.0)
        self.weights.l2 = get_effective_weight(self._recon_l2_cfg, epoch, key="l2_weight", default=0.25)
        self.weights.gradient = get_effective_weight(self._sections["gradient"], epoch, default=0.08)
        self.weights.energy = get_effective_weight(self._sections["energy"], epoch, default=0.02)
        self.weights.spectral = get_effective_weight(self._sections["spectral"], epoch, default=0.02)
        self.weights.radial_spectral = get_effective_weight(self._sections["radial_spectral"], epoch, default=0.02)
        self.weights.subband = get_effective_weight(self._sections["subband"], epoch, default=0.02)
        self.weights.corrupt = get_effective_weight(
            self._sections["latent_corruption"],
            epoch,
            default=0.25,
        )

    def current_weights(self) -> dict[str, float]:
        return {
            "l1": self.weights.l1,
            "l2": self.weights.l2,
            "gradient": self.weights.gradient,
            "energy": self.weights.energy,
            "spectral": self.weights.spectral,
            "radial_spectral": self.weights.radial_spectral,
            "subband": self.weights.subband,
            "corrupt": self.weights.corrupt,
        }

    def _domain_component_scale(self, section_name: str, prefix: str) -> float:
        if not self._scheduled:
            return 1.0
        if prefix.startswith("corrupt_"):
            prefix = prefix.removeprefix("corrupt_")
        domain = "ocn" if prefix.startswith("ocn") else "atm"
        section = self._sections.get(section_name)
        if section is None:
            return 1.0
        return float(section.get(f"{domain}_weight", 1.0))

    def _domain_loss(
        self,
        pred: torch.Tensor,
        target: torch.Tensor,
        prefix: str,
        mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        zero = pred.new_tensor(0.0)
        l1 = _masked_l1_loss(pred, target, mask)
        l2 = _masked_mse_loss(pred, target, mask)
        grad = gradient_loss_2d(pred, target, mask) if self.weights.gradient > 0.0 else zero
        ene = energy_loss(pred, target, mask) if self.weights.energy > 0.0 else zero
        spec, radial, subband = spectral_bundle_loss_2d(
            pred,
            target,
            need_spectral=self.weights.spectral > 0.0,
            need_radial=self.weights.radial_spectral > 0.0,
            need_subband=self.weights.subband > 0.0 and self.subband_mode not in {"spatial", "spatial_band", "palm"},
            bands=self.subband_bands,
            band_weights=self.subband_weights,
            mask=mask,
        )
        if self.weights.subband > 0.0 and self.subband_mode in {"spatial", "spatial_band", "palm"}:
            subband = spatial_subband_loss_2d(
                pred,
                target,
                mask=mask,
                bands=self.subband_bands,
                band_weights=self.subband_weights,
                transition_ratio=self.subband_transition_ratio,
                norm=self.subband_norm,
                relative=self.subband_relative,
                eps=self.subband_eps,
            )
        l1_weight = self.weights.l1
        l2_weight = self.weights.l2
        gradient_weight = self.weights.gradient * self._domain_component_scale("gradient", prefix)
        energy_weight = self.weights.energy * self._domain_component_scale("energy", prefix)
        spectral_weight = self.weights.spectral * self._domain_component_scale("spectral", prefix)
        radial_weight = self.weights.radial_spectral * self._domain_component_scale("radial_spectral", prefix)
        subband_weight = self.weights.subband * self._domain_component_scale("subband", prefix)
        l1_weighted = l1_weight * l1
        l2_weighted = l2_weight * l2
        grad_weighted = gradient_weight * grad
        ene_weighted = energy_weight * ene
        spec_weighted = spectral_weight * spec
        radial_weighted = radial_weight * radial
        subband_weighted = subband_weight * subband
        total = (
            l1_weighted
            + l2_weighted
            + grad_weighted
            + ene_weighted
            + spec_weighted
            + radial_weighted
            + subband_weighted
        )
        return {
            f"{prefix}_l1_loss": l1,
            f"{prefix}_l1_weighted_loss": l1_weighted,
            f"{prefix}_l2_loss": l2,
            f"{prefix}_l2_weighted_loss": l2_weighted,
            f"{prefix}_gradient_loss": grad,
            f"{prefix}_gradient_weighted_loss": grad_weighted,
            f"{prefix}_energy_loss": ene,
            f"{prefix}_energy_weighted_loss": ene_weighted,
            f"{prefix}_spectral_loss": spec,
            f"{prefix}_spectral_weighted_loss": spec_weighted,
            f"{prefix}_radial_spectral_loss": radial,
            f"{prefix}_radial_spectral_weighted_loss": radial_weighted,
            f"{prefix}_subband_loss": subband,
            f"{prefix}_subband_weighted_loss": subband_weighted,
            f"{prefix}_total": total,
        }

    def forward(
        self,
        atm_pred: torch.Tensor,
        ocn_pred: torch.Tensor,
        atm_target: torch.Tensor,
        ocn_target: torch.Tensor,
        atm_corrupt_pred: torch.Tensor | None = None,
        ocn_corrupt_pred: torch.Tensor | None = None,
        atm_mask: torch.Tensor | None = None,
        ocn_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        atm = self._domain_loss(atm_pred, atm_target, "atm", atm_mask)
        ocn = self._domain_loss(ocn_pred, ocn_target, "ocn", ocn_mask)
        total = atm["atm_total"] + ocn["ocn_total"]

        corrupt_total = total.new_tensor(0.0)
        if atm_corrupt_pred is not None and ocn_corrupt_pred is not None:
            corrupt_atm = self._domain_loss(atm_corrupt_pred, atm_target, "corrupt_atm", atm_mask)
            corrupt_ocn = self._domain_loss(ocn_corrupt_pred, ocn_target, "corrupt_ocn", ocn_mask)
            corrupt_total = corrupt_atm["corrupt_atm_total"] + corrupt_ocn["corrupt_ocn_total"]
            total = total + self.weights.corrupt * corrupt_total
        else:
            corrupt_atm = {}
            corrupt_ocn = {}

        losses = {
            **atm,
            **ocn,
            **corrupt_atm,
            **corrupt_ocn,
            "corrupt_total": corrupt_total,
            "total_loss": total,
        }
        return losses
