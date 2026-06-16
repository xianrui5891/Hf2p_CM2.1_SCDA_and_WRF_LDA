from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Union

import numpy as np
import torch
import torch.nn.functional as F
from utils.plotting import FieldPlotter

from model import CoupledAE
from utils.data import RUNTIME_ATM_VARS, RUNTIME_OCN_VARS
from utils.normalization import NormalizationConfig, VariableNormalizer
from utils.paths import (
    DEFAULT_METADATA_PATH,
    DEFAULT_NMC_BACKGROUND_COVARIANCE_PATH,
    DEFAULT_OBS_PATH,
    DEFAULT_AE_MODEL_PATH,
    PROJECT_DIR,
)


class latent_space_da:
    """DA adapter kept compatible with temp/cm2_cda_main.py."""

    def __init__(
        self,
        da_interval: int,
        dev: str | torch.device = "cuda",
        learning_rate: float = 1e-3,
        max_iterations: int = 30,
        atm_lat_size: int = 90,
        atm_lon_size: int = 144,
        ocn_lat_size: int = 200,
        ocn_lon_size: int = 360,
        obs_path: str | Path = DEFAULT_OBS_PATH,
        ae_model_path: str | Path = DEFAULT_AE_MODEL_PATH,
        metadata_path: str | Path = DEFAULT_METADATA_PATH,
        nmc_background_covariance_path: str | Path | None = DEFAULT_NMC_BACKGROUND_COVARIANCE_PATH,
        b_diag_floor: float = 1e-12,
        ps_obs_error_variance: float = 1e-2,
        sst_obs_error_variance: float = 1e-2,
        normalize_3dvar_cost: bool = True,
        coastal_ocean_relaxation: bool = False,
        coastal_ocean_relaxation_width: int = 2,
        coastal_ocean_relaxation_factor: float = 0.35,
        coastal_sst_obs_weighting: bool = False,
        coastal_sst_obs_weight_width: int = 2,
        coastal_sst_obs_weight: float = 2.0,
        da_diagnostic_plot_dir: str | Path | None = "python_plots/diagnostics",
    ) -> None:
        del atm_lat_size, atm_lon_size, ocn_lat_size, ocn_lon_size
        self.interval = int(da_interval)
        self.dev = torch.device(dev)
        self.learning_rate = float(learning_rate)
        self.max_iterations = int(max_iterations)
        self.ps_obs_error_variance = max(float(ps_obs_error_variance), 1e-12)
        self.sst_obs_error_variance = max(float(sst_obs_error_variance), 1e-12)
        self.normalize_3dvar_cost = bool(normalize_3dvar_cost)
        self.coastal_ocean_relaxation = bool(coastal_ocean_relaxation)
        self.coastal_ocean_relaxation_width = max(int(coastal_ocean_relaxation_width), 0)
        self.coastal_ocean_relaxation_factor = float(np.clip(coastal_ocean_relaxation_factor, 0.0, 1.0))
        self.coastal_sst_obs_weighting = bool(coastal_sst_obs_weighting)
        self.coastal_sst_obs_weight_width = max(int(coastal_sst_obs_weight_width), 0)
        self.coastal_sst_obs_weight = max(float(coastal_sst_obs_weight), 1.0)
        self.da_diagnostic_plotter = FieldPlotter(da_diagnostic_plot_dir) if da_diagnostic_plot_dir else None

        obs_path = obs_path or DEFAULT_OBS_PATH
        ae_model_path = ae_model_path or DEFAULT_AE_MODEL_PATH
        metadata_path = metadata_path or DEFAULT_METADATA_PATH
        nmc_background_covariance_path = nmc_background_covariance_path or DEFAULT_NMC_BACKGROUND_COVARIANCE_PATH

        self.obs = torch.load(obs_path, map_location="cpu", weights_only=False)
        self.obs_data = self.obs["data"]
        with open(metadata_path, "r", encoding="utf-8") as f:
            self.metadata = json.load(f)
        self.normalizer = VariableNormalizer(NormalizationConfig(metadata_path=metadata_path))
        self.atm_vars = list(self.metadata.get("atm_vars", RUNTIME_ATM_VARS))
        self.ocn_vars = list(self.metadata.get("ocn_vars", RUNTIME_OCN_VARS))

        raw_ps = self.obs_data["ps1"].squeeze().to(self.dev)
        raw_sst = self.obs_data["sst"].squeeze().to(self.dev)
        self.obs_masks = {
            "ps1": torch.isfinite(raw_ps),
            "sst": torch.isfinite(raw_sst),
        }
        self.obs_data["ps1"] = self.normalizer.normalize(
            raw_ps,
            var_name="PS",
            domain="atm",
            method="standardize",
        ).masked_fill(~self.obs_masks["ps1"], float("nan"))
        self.obs_data["sst"] = self.normalizer.normalize(
            raw_sst,
            var_name="SST",
            domain="ocn",
            method="standardize",
        ).masked_fill(~self.obs_masks["sst"], float("nan"))

        self.ae = self._load_model(ae_model_path).to(self.dev).eval()
        self.latent_dim = int(getattr(self.ae, "latent_dim", 0))
        self.B_diag = self._load_background_covariance_diag(
            nmc_background_covariance_path,
            expected_dim=self.latent_dim,
            b_diag_floor=b_diag_floor,
        )

    def _load_model(self, checkpoint_path: str | Path) -> CoupledAE:
        ckpt = torch.load(checkpoint_path, map_location=self.dev, weights_only=False)
        model_cfg = ckpt.get("model_config", {}) if isinstance(ckpt, dict) else {}
        latent_dim = ckpt.get("latent_dim", model_cfg.get("latent_dim")) if isinstance(ckpt, dict) else None
        state = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt)) if isinstance(ckpt, dict) else ckpt
        state = {k.removeprefix("module."): v for k, v in state.items()}
        fusion_mode = model_cfg.get("fusion_mode")
        if not fusion_mode:
            coupling_keys = [key for key in state if ".pre_down_coupling." in key or ".post_down_coupling." in key]
            fusion_mode = "cnn" if any(".net." in key for key in coupling_keys) else "cross_attention"
        width = int(model_cfg.get("width", 64))
        model = CoupledAE(
            latent_dim=latent_dim,
            latent_channels=model_cfg.get("latent_channels"),
            in_atm_channels=model_cfg.get("in_atm_channels", 4),
            in_ocn_channels=model_cfg.get("in_ocn_channels", 4),
            width=width,
            decoder_width=model_cfg.get("decoder_width", width * 2),
            decoder_middle_depth=model_cfg.get("decoder_middle_depth", 2),
            decoder_stage_depth=model_cfg.get("decoder_stage_depth", 1),
            atm_decoder_stage_depth=model_cfg.get("atm_decoder_stage_depth"),
            ocn_decoder_stage_depth=model_cfg.get("ocn_decoder_stage_depth"),
            cross_depth=model_cfg.get("cross_depth", 4),
            heads=model_cfg.get("heads", 4),
            attention_backend=model_cfg.get("attention_backend", "auto"),
            fusion_mode=fusion_mode,
            cnn_fusion_depth=model_cfg.get("cnn_fusion_depth", 2),
        )
        missing, unexpected = model.load_state_dict(state, strict=False)
        if missing:
            print(f"[WARN] Missing model keys: {len(missing)}")
        if unexpected:
            print(f"[WARN] Unexpected model keys: {len(unexpected)}")
        print(f"AE model loaded from {checkpoint_path}. latent_dim={model.latent_dim}")
        return model

    def _resolve_path(self, path: Optional[Union[str, Path]]) -> Optional[Path]:
        if path is None:
            return None
        path = Path(path).expanduser()
        if not path.is_absolute():
            path = PROJECT_DIR / path
        return path.resolve()

    def _load_background_covariance_diag(
        self,
        covariance_path: Optional[Union[str, Path]],
        expected_dim: int,
        b_diag_floor: float,
    ) -> torch.Tensor:
        fallback = torch.ones(expected_dim, dtype=torch.float32, device=self.dev)
        path = self._resolve_path(covariance_path)
        if path is None or not path.exists():
            print("[WARN] NMC covariance is unavailable. Using unit B_diag.")
            return fallback
        with np.load(path, allow_pickle=True) as data:
            key = next((item for item in ("B_diag", "background_covariance_diag", "variance") if item in data), None)
            if key is None:
                raise KeyError(f"{path} must contain B_diag/background_covariance_diag/variance.")
            diag = np.asarray(data[key], dtype=np.float32).reshape(-1)
        if diag.size != expected_dim:
            raise ValueError(f"NMC B_diag length mismatch: got {diag.size}, expected {expected_dim}.")
        floor = max(float(b_diag_floor), 1e-12)
        diag = np.clip(np.nan_to_num(diag, nan=floor, posinf=np.finfo(np.float32).max, neginf=floor), floor, None)
        print(f"[INFO] Loaded NMC B_diag from {path}.")
        return torch.as_tensor(diag, dtype=torch.float32, device=self.dev)

    @torch.no_grad()
    def to_tensor_from_list(self, type: str, np_list, device: str | torch.device | None = None) -> torch.Tensor:
        vars_ = self.atm_vars if type == "atm" else self.ocn_vars
        if len(np_list) != len(vars_):
            raise ValueError(f"{type} channel count {len(np_list)} does not match vars {len(vars_)}.")
        normalized = []
        for arr, var_name in zip(np_list, vars_):
            item = np.ma.filled(arr, np.nan) if np.ma.isMaskedArray(arr) else arr
            item = self.normalizer.normalize(item, var_name=var_name, domain=type, method="standardize")
            if isinstance(item, torch.Tensor):
                item = item.detach().cpu().numpy()
            normalized.append(np.asarray(item, dtype=np.float32).T)
        tensor = torch.from_numpy(np.stack(normalized, axis=0)).unsqueeze(0)
        return tensor.to(device or self.dev)

    def mask_from_list(self, type: str, np_list, device: str | torch.device | None = None) -> torch.Tensor:
        vars_ = self.atm_vars if type == "atm" else self.ocn_vars
        if len(np_list) != len(vars_):
            raise ValueError(f"{type} channel count {len(np_list)} does not match vars {len(vars_)}.")
        masks = []
        for arr in np_list:
            if np.ma.isMaskedArray(arr):
                valid = ~np.ma.getmaskarray(arr)
                values = np.ma.filled(arr, np.nan)
                valid = valid & np.isfinite(np.asarray(values, dtype=np.float32))
            else:
                valid = np.isfinite(np.asarray(arr, dtype=np.float32))
            masks.append(np.asarray(valid, dtype=bool).T)
        return torch.from_numpy(np.stack(masks, axis=0)).unsqueeze(0).to(device or self.dev)

    def _broadcast_surface_mask(self, mask2d: torch.Tensor, channels: int) -> torch.Tensor:
        mask2d = mask2d.to(self.dev, dtype=torch.bool)
        return mask2d.unsqueeze(0).unsqueeze(0).expand(1, channels, *mask2d.shape)

    def _output_valid_masks(
        self,
        nc_idx: int,
        mask_atm_input: torch.Tensor,
        mask_ocn_input: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        atm_valid = mask_atm_input.to(self.dev, dtype=torch.bool)
        ocn_valid = mask_ocn_input.to(self.dev, dtype=torch.bool)
        sst_valid = self.obs_masks["sst"][nc_idx] & torch.isfinite(self.obs_data["sst"][nc_idx])
        ocn_valid = ocn_valid & self._broadcast_surface_mask(sst_valid, ocn_valid.shape[1])
        return atm_valid, ocn_valid

    def _coastal_surface_mask(self, ocean: torch.Tensor, width: int) -> torch.Tensor:
        ocean = ocean.to(self.dev, dtype=torch.bool)
        if ocean.dim() != 2 or width <= 0:
            return torch.zeros_like(ocean, dtype=torch.bool)
        near_land = (~ocean).unsqueeze(0).unsqueeze(0).to(dtype=torch.float32)
        for _ in range(width):
            near_land = F.max_pool2d(near_land, kernel_size=3, stride=1, padding=1)
        return ((near_land.squeeze(0).squeeze(0) > 0.0) & ocean)

    def _coastal_ocean_increment_weight(self, output_ocn_valid: torch.Tensor) -> torch.Tensor | None:
        if (
            not self.coastal_ocean_relaxation
            or self.coastal_ocean_relaxation_width <= 0
            or self.coastal_ocean_relaxation_factor >= 1.0
        ):
            return None
        sst_idx = self.ocn_vars.index("SST") if "SST" in self.ocn_vars else 0
        ocean = output_ocn_valid[:, sst_idx : sst_idx + 1].to(self.dev, dtype=torch.bool)
        coastal = self._coastal_surface_mask(ocean[0, 0], self.coastal_ocean_relaxation_width).unsqueeze(0).unsqueeze(0)
        coastal = coastal & ocean
        weight = torch.ones_like(ocean, dtype=torch.float32)
        weight = torch.where(
            coastal,
            torch.full_like(weight, self.coastal_ocean_relaxation_factor),
            weight,
        )
        return weight.expand_as(output_ocn_valid).to(device=self.dev)

    def _relax_coastal_ocean_increment(
        self,
        input_ocn_t: torch.Tensor,
        ocn_tensor: torch.Tensor,
        output_ocn_valid: torch.Tensor,
    ) -> torch.Tensor:
        weight = self._coastal_ocean_increment_weight(output_ocn_valid)
        if weight is None:
            return ocn_tensor
        coastal_points = int(((weight[:, :1] < 1.0) & output_ocn_valid[:, :1]).sum().item())
        print(
            f"[DA] coastal ocean relaxation width={self.coastal_ocean_relaxation_width} "
            f"factor={self.coastal_ocean_relaxation_factor:g} points={coastal_points}"
        )
        relaxed = input_ocn_t + (ocn_tensor - input_ocn_t) * weight.to(dtype=ocn_tensor.dtype)
        return torch.where(output_ocn_valid, relaxed, ocn_tensor)

    @torch.no_grad()
    def _plot_sst_da_diagnostics(
        self,
        nc: int,
        nc_idx: int,
        mask_ocn: torch.Tensor,
        coastal_sst: torch.Tensor | None,
        input_ocn_t: torch.Tensor,
        background_latent: torch.Tensor,
        analysis_latent: torch.Tensor,
    ) -> None:
        if self.da_diagnostic_plotter is None:
            return
        sst_idx = self.ocn_vars.index("SST") if "SST" in self.ocn_vars else 0
        _, decoded_bg_ocn = self.ae.decode(background_latent)
        _, an_ocn = self.ae.decode(analysis_latent)
        obs_norm = self.obs_data["sst"][nc_idx]
        input_bg_norm = input_ocn_t[0, sst_idx]
        decoded_bg_norm = decoded_bg_ocn[0, sst_idx]
        an_norm = an_ocn[0, sst_idx]

        def _denorm(field: torch.Tensor) -> np.ndarray:
            arr = self.normalizer.denormalize(field.detach().cpu(), var_name="SST", domain="ocn", method="standardize")
            if isinstance(arr, torch.Tensor):
                arr = arr.numpy()
            return np.asarray(arr, dtype=np.float64).T

        valid = mask_ocn.detach().cpu().numpy().T
        obs = _denorm(obs_norm)
        input_bg = _denorm(input_bg_norm)
        decoded_bg = _denorm(decoded_bg_norm)
        an = _denorm(an_norm)
        obs_plot = np.ma.array(obs, mask=~valid)
        input_bg_residual = np.ma.array(obs - input_bg, mask=~valid)
        decoded_bg_residual = np.ma.array(obs - decoded_bg, mask=~valid)
        an_residual = np.ma.array(obs - an, mask=~valid)
        improvement = np.ma.array(np.abs(input_bg_residual) - np.abs(an_residual), mask=~valid)
        time_str = f"{nc}_obs{nc_idx:03d}"
        self.da_diagnostic_plotter.plot_and_save(
            time_str=time_str,
            var_name="sst_observation",
            data2d=obs_plot,
            data_type="diagnostics",
        )
        self.da_diagnostic_plotter.plot_and_save(
            time_str=time_str,
            var_name="sst_background_residual",
            data2d=input_bg_residual,
            data_type="diagnostics",
        )
        self.da_diagnostic_plotter.plot_and_save(
            time_str=time_str,
            var_name="sst_decoded_background_residual",
            data2d=decoded_bg_residual,
            data_type="diagnostics",
        )
        self.da_diagnostic_plotter.plot_and_save(
            time_str=time_str,
            var_name="sst_analysis_residual",
            data2d=an_residual,
            data_type="diagnostics",
        )
        self.da_diagnostic_plotter.plot_and_save(
            time_str=time_str,
            var_name="sst_abs_residual_improvement",
            data2d=improvement,
            data_type="diagnostics",
        )
        if coastal_sst is not None:
            coastal = coastal_sst.detach().cpu().numpy().T.astype(np.float64)
            coastal = np.ma.array(coastal, mask=~valid)
            self.da_diagnostic_plotter.plot_and_save(
                time_str=time_str,
                var_name="sst_coastal_obs_mask",
                data2d=coastal,
                data_type="diagnostics",
                vmin=0.0,
                vmax=1.0,
            )

    @torch.no_grad()
    def tensor_to_numpy_list(
        self,
        type: str,
        tensor: torch.Tensor,
        valid_mask: torch.Tensor | None = None,
        templates: list | None = None,
    ) -> list[np.ndarray]:
        vars_ = self.atm_vars if type == "atm" else self.ocn_vars
        tensor = tensor.detach().cpu()
        if tensor.dim() == 4 and tensor.shape[0] == 1:
            tensor = tensor.squeeze(0)
        if valid_mask is not None:
            valid_mask = valid_mask.detach().cpu()
            if valid_mask.dim() == 4 and valid_mask.shape[0] == 1:
                valid_mask = valid_mask.squeeze(0)
        out = []
        for idx, var_name in enumerate(vars_):
            arr = self.normalizer.denormalize(tensor[idx], var_name=var_name, domain=type, method="standardize")
            if isinstance(arr, torch.Tensor):
                arr = arr.numpy()
            arr_np = np.asarray(arr, dtype=np.float64).T
            template = templates[idx] if templates is not None and idx < len(templates) else None
            valid = None
            if valid_mask is not None:
                valid = np.asarray(valid_mask[idx].numpy(), dtype=bool).T
                if template is not None:
                    template_values = np.ma.filled(template, np.nan) if np.ma.isMaskedArray(template) else template
                    template_values = np.asarray(template_values, dtype=np.float64)
                    if template_values.shape == arr_np.shape:
                        arr_np = np.where(valid, arr_np, template_values)
                    else:
                        arr_np = np.where(valid, arr_np, np.nan)
                else:
                    arr_np = np.where(valid, arr_np, np.nan)
            if np.ma.isMaskedArray(template):
                template_mask = np.ma.getmaskarray(template)
                if valid is not None:
                    template_mask = template_mask | ~valid
                out.append(np.ma.array(arr_np, mask=template_mask, copy=False))
            else:
                out.append(arr_np)
        return out

    def mask_field_for_plot(self, var_name: str, data2d, nc: int | None = None):
        name = var_name.lower()
        if name in {"sst", "sea_surface_temperature"}:
            idx = int(nc / self.interval) if nc is not None else 0
            idx = max(0, min(idx, int(self.obs_masks["sst"].shape[0]) - 1))
            valid = (self.obs_masks["sst"][idx] & torch.isfinite(self.obs_data["sst"][idx])).detach().cpu().numpy().T
            arr = np.asarray(np.ma.filled(data2d, np.nan) if np.ma.isMaskedArray(data2d) else data2d)
            if arr.shape == valid.shape:
                return np.ma.array(arr, mask=~valid, copy=False)
        if name in {"p_bot", "ps", "ps1"}:
            idx = int(nc / self.interval) if nc is not None else 0
            idx = max(0, min(idx, int(self.obs_masks["ps1"].shape[0]) - 1))
            valid = (self.obs_masks["ps1"][idx] & torch.isfinite(self.obs_data["ps1"][idx])).detach().cpu().numpy().T
            arr = np.asarray(np.ma.filled(data2d, np.nan) if np.ma.isMaskedArray(data2d) else data2d)
            if arr.shape == valid.shape:
                return np.ma.array(arr, mask=~valid, copy=False)
        return data2d

    def H(self, latent_x: torch.Tensor, obs_mask: list[torch.Tensor]) -> torch.Tensor:
        obs_atm, obs_ocn = self.ae.decode(latent_x)
        obs_atm = torch.nan_to_num(obs_atm.squeeze(0), nan=0.0, posinf=0.0, neginf=0.0)
        obs_ocn = torch.nan_to_num(obs_ocn.squeeze(0), nan=0.0, posinf=0.0, neginf=0.0)
        mask_atm = obs_mask[0].to(obs_atm.device)
        mask_ocn = obs_mask[1].to(obs_ocn.device)
        ps_idx = self.atm_vars.index("PS") if "PS" in self.atm_vars else 3
        sst_idx = self.ocn_vars.index("SST") if "SST" in self.ocn_vars else 0
        return torch.cat((obs_atm[ps_idx][mask_atm], obs_ocn[sst_idx][mask_ocn]), dim=0).squeeze()

    @torch.enable_grad()
    def apply_3DVar(
        self,
        H,
        B_diag: torch.Tensor,
        R_diag: torch.Tensor,
        xb: torch.Tensor,
        y: torch.Tensor,
        obs_mask: list[torch.Tensor],
        max_iterations: int = 100,
        learning_rate: float = 1e-3,
    ) -> torch.Tensor:
        x = torch.nn.Parameter(xb.detach().clone())
        opt = torch.optim.Adam([x], lr=learning_rate)
        B_diag = B_diag.to(x.device)
        R_diag = R_diag.to(x.device)
        best_x = x.detach().clone()
        for _ in range(max_iterations):
            opt.zero_grad(set_to_none=True)
            dx = x.ravel() - xb.ravel()
            background_terms = dx.square() / B_diag
            loss_b = background_terms.mean() if self.normalize_3dvar_cost else background_terms.sum()
            hx = H(x, obs_mask).ravel()
            residual = y.ravel() - hx
            valid_residual = torch.isfinite(residual) & torch.isfinite(R_diag) & (R_diag > 0)
            if residual.numel() != R_diag.numel():
                raise ValueError(f"H(x) size {residual.numel()} does not match R_diag size {R_diag.numel()}.")
            if valid_residual.any():
                residual = residual[valid_residual]
                R_valid = R_diag[valid_residual]
                obs_terms = residual.square() / R_valid
                loss_o = obs_terms.mean() if self.normalize_3dvar_cost else obs_terms.sum()
            else:
                loss_o = torch.zeros((), dtype=x.dtype, device=x.device)
            loss = loss_b + loss_o
            loss.backward()
            opt.step()
            best_x = x.detach().clone()
        return best_x

    @torch.no_grad()
    def do_da(self, nc: int, input_atm: list, input_ocn: list):
        mask_atm_input = self.mask_from_list("atm", input_atm)
        mask_ocn_input = self.mask_from_list("ocn", input_ocn)
        input_atm_t = self.to_tensor_from_list("atm", input_atm)
        input_ocn_t = self.to_tensor_from_list("ocn", input_ocn)
        input_atm_t = torch.nan_to_num(input_atm_t, nan=0.0, posinf=0.0, neginf=0.0)
        input_ocn_t = torch.nan_to_num(input_ocn_t, nan=0.0, posinf=0.0, neginf=0.0)

        latent_x = self.ae.encode(input_atm_t, input_ocn_t)
        nc_idx = int(nc / self.interval)
        ps_idx = self.atm_vars.index("PS") if "PS" in self.atm_vars else 3
        sst_idx = self.ocn_vars.index("SST") if "SST" in self.ocn_vars else 0
        mask_atm = (
            self.obs_masks["ps1"][nc_idx]
            & torch.isfinite(self.obs_data["ps1"][nc_idx])
            & mask_atm_input[0, ps_idx]
        ).to(self.dev)
        mask_ocn = (
            self.obs_masks["sst"][nc_idx]
            & torch.isfinite(self.obs_data["sst"][nc_idx])
            & mask_ocn_input[0, sst_idx]
        ).to(self.dev)
        obs_values = torch.cat(
            (
                self.obs_data["ps1"][nc_idx][mask_atm],
                self.obs_data["sst"][nc_idx][mask_ocn],
            ),
            dim=0,
        ).to(self.dev)
        output_atm_valid, output_ocn_valid = self._output_valid_masks(nc_idx, mask_atm_input, mask_ocn_input)
        if obs_values.numel() == 0:
            atm_tensor = input_atm_t.masked_fill(~output_atm_valid, float("nan"))
            ocn_tensor = input_ocn_t.masked_fill(~output_ocn_valid, float("nan"))
            return (
                self.tensor_to_numpy_list("atm", atm_tensor, valid_mask=output_atm_valid, templates=input_atm),
                self.tensor_to_numpy_list("ocn", ocn_tensor, valid_mask=output_ocn_valid, templates=input_ocn),
            )
        ps_count = int(mask_atm.sum().item())
        sst_count = int(mask_ocn.sum().item())
        sst_r_diag = torch.full((sst_count,), self.sst_obs_error_variance, dtype=torch.float32, device=self.dev)
        coastal_sst_obs_mask = torch.zeros((sst_count,), dtype=torch.bool, device=self.dev)
        coastal_sst = None
        coastal_sst_count = 0
        if self.coastal_sst_obs_weighting and self.coastal_sst_obs_weight_width > 0 and sst_count > 0:
            coastal_sst = self._coastal_surface_mask(mask_ocn, self.coastal_sst_obs_weight_width)
            coastal_sst_obs_mask = coastal_sst[mask_ocn]
            coastal_sst_count = int((coastal_sst & mask_ocn).sum().item())
            if coastal_sst_count:
                sst_r_diag = torch.where(
                    coastal_sst_obs_mask,
                    sst_r_diag / self.coastal_sst_obs_weight,
                    sst_r_diag,
                )
        R_diag = torch.cat(
            (
                torch.full((ps_count,), self.ps_obs_error_variance, dtype=torch.float32, device=self.dev),
                sst_r_diag,
            ),
            dim=0,
        )
        assim_latent_x = self.apply_3DVar(
            self.H,
            self.B_diag,
            R_diag,
            latent_x,
            obs_values,    
            [mask_atm, mask_ocn],
            self.max_iterations,
            self.learning_rate,
        )
        self._plot_sst_da_diagnostics(nc, nc_idx, mask_ocn, coastal_sst, input_ocn_t, latent_x, assim_latent_x)
        atm_tensor, ocn_tensor = self.ae.decode(assim_latent_x)
        ocn_tensor = self._relax_coastal_ocean_increment(input_ocn_t, ocn_tensor, output_ocn_valid)
        atm_tensor = atm_tensor.masked_fill(~output_atm_valid, float("nan"))
        ocn_tensor = ocn_tensor.masked_fill(~output_ocn_valid, float("nan"))
        return (
            self.tensor_to_numpy_list("atm", atm_tensor, valid_mask=output_atm_valid, templates=input_atm),
            self.tensor_to_numpy_list("ocn", ocn_tensor, valid_mask=output_ocn_valid, templates=input_ocn),
        )
