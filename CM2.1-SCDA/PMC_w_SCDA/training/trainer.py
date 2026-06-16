from __future__ import annotations

import gc
import json
import math
import time
from contextlib import nullcontext
from pathlib import Path
from typing import Any

import torch
from torch.amp import GradScaler, autocast
from torch.nn.parallel import DistributedDataParallel
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import DataLoader

from model import CoupledAE
from training.checkpoint import save_checkpoint, unwrap_model
from training.losses import CoupledAELoss, get_loss_section, get_schedule_factor
from training.optimizers import build_optimizer, build_scheduler
from utils.ema import EMAModel
from utils.logger import get_logger
from utils.plotting import LossPlotter


class CoupledAETrainer:
    """PALM-GDA-style trainer infra specialized only at the model/loss boundary."""

    def __init__(
        self,
        cfg: dict[str, Any],
        train_loader: DataLoader,
        val_loader: DataLoader,
        dataset: Any,
        device: torch.device,
        backend: str = "none",
        distributed: bool = False,
        local_rank: int = 0,
    ) -> None:
        self.cfg = cfg
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.dataset = dataset
        self.device = device
        self.backend = str(backend).lower()
        self.distributed = distributed
        self.local_rank = local_rank
        self.logger = get_logger()
        self.ds_engine = None

        model_cfg = cfg.get("model", {})
        self.model = CoupledAE(
            latent_channels=int(model_cfg.get("latent_channels", 36)),
            width=int(model_cfg.get("width", 56)),
            decoder_width=int(model_cfg.get("decoder_width", 160)),
            decoder_middle_depth=int(model_cfg.get("decoder_middle_depth", 4)),
            decoder_stage_depth=model_cfg.get("decoder_stage_depth", 2),
            atm_decoder_stage_depth=model_cfg.get("atm_decoder_stage_depth"),
            ocn_decoder_stage_depth=model_cfg.get("ocn_decoder_stage_depth"),
            cross_depth=int(model_cfg.get("cross_depth", 4)),
            heads=int(model_cfg.get("heads", 4)),
            attention_backend=cfg.get("backends", {}).get("attention", model_cfg.get("attention_backend", "auto")),
            fusion_mode=model_cfg.get("fusion_mode", "cross_attention"),
            cnn_fusion_depth=int(model_cfg.get("cnn_fusion_depth", 2)),
        ).to(device)
        encoder_params = sum(p.numel() for p in self.model.encoder.parameters())
        decoder_params = sum(p.numel() for p in self.model.atm_decoder.parameters()) + sum(p.numel() for p in self.model.ocn_decoder.parameters())
        total_params = sum(p.numel() for p in self.model.parameters())
        self.logger.info(
            "Model initialized: latent_shape=%s latent_dim=%d total_params=%d encoder_params=%d decoder_params=%d decoder/encoder=%.3f",
            self.model.latent_shape,
            self.model.latent_dim,
            total_params,
            encoder_params,
            decoder_params,
            decoder_params / max(encoder_params, 1),
        )
        self._apply_runtime_optimizations()
        self.ddp_model = None
        if self.backend == "ddp" and distributed:
            self.model = DistributedDataParallel(self.model, device_ids=[local_rank] if torch.cuda.is_available() else None)
            self.ddp_model = self.model

        loss_cfg = cfg.get("loss", {})
        self.criterion = CoupledAELoss(loss_cfg=loss_cfg)
        self.latent_corruption_cfg = get_loss_section(loss_cfg, "latent_corruption")
        if "weight" not in self.latent_corruption_cfg and "corrupt" in loss_cfg:
            self.latent_corruption_cfg["weight"] = loss_cfg["corrupt"]
        self.corruption_enabled = bool(self.latent_corruption_cfg.get("enabled", True))
        self.corruption_mode = str(self.latent_corruption_cfg.get("mode", "single_pass_denoising")).lower()
        self.base_corruption_sigma = float(
            self.latent_corruption_cfg.get("sigma", cfg.get("training", {}).get("latent_corruption_sigma", 0.03))
        )
        self.current_corruption_sigma = self.base_corruption_sigma if self.corruption_enabled else 0.0
        logging_cfg = cfg.get("logging", {})
        self.slow_data_threshold_s = float(logging_cfg.get("slow_data_threshold_s", 30.0))
        self.slow_compute_threshold_s = float(logging_cfg.get("slow_compute_threshold_s", 30.0))

        train_cfg = cfg.get("training", {})
        self.epochs = int(train_cfg.get("epochs", 20))
        self.amp_enabled = bool(train_cfg.get("amp_enabled", torch.cuda.is_available()))
        amp_dtype = str(train_cfg.get("amp_dtype", "float16")).lower()
        self.amp_dtype = torch.bfloat16 if "bf16" in amp_dtype or "bfloat16" in amp_dtype else torch.float16
        self.scaler = GradScaler("cuda", enabled=self.amp_enabled and self.amp_dtype == torch.float16)
        self.gradient_accumulation_steps = max(int(train_cfg.get("gradient_accumulation_steps", 1)), 1)
        self.gradient_clip_norm = float(train_cfg.get("gradient_clip_norm", 1.0))
        self.log_every_n_steps = int(self.cfg.get("logging", {}).get("log_every_n_steps", 10))
        self.val_every_n_epochs = max(
            int(self.cfg.get("logging", {}).get("val_every_n_epochs", train_cfg.get("val_every_n_epochs", 1))),
            1,
        )
        opt_cfg = dict(train_cfg.get("optimizer", {}))
        opt_cfg.setdefault("name", "adamw")
        opt_cfg.setdefault("lr", train_cfg.get("lr", 2e-4))
        opt_cfg.setdefault("weight_decay", train_cfg.get("weight_decay", 1e-4))
        self.optimizer = build_optimizer(self.model, opt_cfg)
        sched_cfg = dict(train_cfg.get("scheduler", {}))
        if not sched_cfg:
            sched_cfg = {
                "name": "plateau",
                "factor": train_cfg.get("scheduler_factor", 0.5),
                "patience": train_cfg.get("scheduler_patience", 5),
            }
        self.scheduler_warmup_epochs = int(sched_cfg.get("warmup_epochs", 0))
        self.scheduler_warmup_start_factor = float(sched_cfg.get("warmup_start_factor", 0.1))
        self._base_lrs = [float(group["lr"]) for group in self.optimizer.param_groups]
        self._active_scheduler_epoch = 0
        self.scheduler = build_scheduler(
            self.optimizer,
            sched_cfg,
            total_epochs=max(self.epochs - self.scheduler_warmup_epochs, 1),
            steps_per_epoch=len(train_loader),
        )
        self._apply_warmup_lr(0)
        if self.backend == "deepspeed" and distributed:
            self._init_deepspeed()
        ema_cfg = train_cfg.get("ema", {})
        self.ema: EMAModel | None = None
        if bool(ema_cfg.get("enabled", True)):
            self.ema = EMAModel(
                unwrap_model(self.model),
                decay=float(ema_cfg.get("decay", 0.9999)),
                warmup_steps=int(ema_cfg.get("warmup_steps", 1000)),
            )
        es_cfg = train_cfg.get("early_stopping", {})
        self.early_stopping_enabled = bool(es_cfg.get("enabled", train_cfg.get("early_stopping_enabled", True)))
        self.early_patience = int(es_cfg.get("patience", train_cfg.get("early_stopping_patience", 30)))
        self.early_min_delta = float(es_cfg.get("min_delta", train_cfg.get("early_stopping_min_delta", 1e-6)))

        self.output_dir = Path(cfg.get("output_dir", "runs/run_ae"))
        ckpt_cfg = cfg.get("checkpoint", {})
        self.checkpoint_dir = Path(ckpt_cfg.get("save_dir", self.output_dir))
        self.save_every_n_epochs = max(int(ckpt_cfg.get("save_every_n_epochs", 1)), 1)
        self.keep_last_n = int(ckpt_cfg.get("keep_last_n", 50))
        self.log_dir = Path(cfg.get("logging", {}).get("log_dir", self.output_dir / "logs"))
        self.plotter = LossPlotter(self.log_dir)
        self.history: list[dict[str, float]] = []
        self.best_val = float("inf")
        self.bad_epochs = 0
        self.start_epoch = 1
        self.global_step = 0
        self._pending_resume_path: Path | None = None
        self._pending_resume_ckpt: dict[str, Any] | None = None
        self.logger.info(
            "Trainer ready: epochs=%d optimizer=%s fused=%s lr=%g batch_size=%s grad_accum=%d val_every=%d amp=%s amp_dtype=%s corruption_mode=%s corruption_sigma=%g checkpoint_dir=%s loss_plot_dir=%s",
            self.epochs,
            opt_cfg.get("name", "adamw"),
            bool(opt_cfg.get("fused", False)),
            float(opt_cfg.get("lr", 2e-4)),
            train_cfg.get("batch_size", 2),
            self.gradient_accumulation_steps,
            self.val_every_n_epochs,
            self.amp_enabled,
            self.amp_dtype,
            self.corruption_mode,
            self.current_corruption_sigma,
            self.checkpoint_dir,
            self.log_dir,
        )

    def _apply_runtime_optimizations(self) -> None:
        runtime_cfg = self.cfg.get("runtime", {})
        backend_cfg = self.cfg.get("backends", {})
        if torch.cuda.is_available():
            torch.backends.cuda.matmul.allow_tf32 = bool(runtime_cfg.get("allow_tf32", False))
            torch.backends.cudnn.allow_tf32 = bool(runtime_cfg.get("allow_tf32", False))
        torch.backends.cudnn.benchmark = bool(runtime_cfg.get("cudnn_benchmark", False))
        precision = runtime_cfg.get("float32_matmul_precision")
        if precision:
            torch.set_float32_matmul_precision(str(precision))
        conv_backend = str(backend_cfg.get("conv", runtime_cfg.get("conv", "torch"))).lower()
        if conv_backend in {"channels_last", "super_conv", "super-conv"} or runtime_cfg.get("channels_last", False):
            self.model = self.model.to(memory_format=torch.channels_last)
            self.channels_last = True
        else:
            self.channels_last = False
        if bool(runtime_cfg.get("torch_compile", False)):
            mode = str(runtime_cfg.get("torch_compile_mode", "reduce-overhead"))
            self.model = torch.compile(self.model, mode=mode)
        if runtime_cfg.get("log_backend_status", True) and self.local_rank == 0:
            self.logger.info(
                "Runtime backends: attention=%s conv=%s channels_last=%s torch_compile=%s tf32=%s",
                backend_cfg.get("attention", "auto"),
                conv_backend,
                self.channels_last,
                bool(runtime_cfg.get("torch_compile", False)),
                bool(runtime_cfg.get("allow_tf32", False)),
            )

    def _log_cuda_memory(self, tag: str, epoch: int) -> None:
        if self.local_rank != 0 or self.device.type != "cuda" or not torch.cuda.is_available():
            return
        allocated = torch.cuda.memory_allocated(self.device) / 1024**3
        reserved = torch.cuda.memory_reserved(self.device) / 1024**3
        peak = torch.cuda.max_memory_allocated(self.device) / 1024**3
        self.logger.info(
            "CUDA memory %s epoch=%03d allocated=%.2fGiB reserved=%.2fGiB peak=%.2fGiB",
            tag,
            epoch,
            allocated,
            reserved,
            peak,
        )

    def _release_epoch_memory(self, tag: str, epoch: int) -> None:
        self._log_cuda_memory(f"{tag}-before-empty-cache", epoch)
        gc.collect()
        if self.device.type == "cuda" and torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(self.device)
        self._log_cuda_memory(f"{tag}-after-empty-cache", epoch)

    def _init_deepspeed(self) -> None:
        try:
            import deepspeed
        except ImportError as exc:
            raise RuntimeError("DeepSpeed backend requested but deepspeed is not installed.") from exc

        ds_cfg = self.cfg.get("distributed", {}).get("deepspeed", {})
        train_cfg = self.cfg.get("training", {})
        micro_batch = int(ds_cfg.get("train_micro_batch_size_per_gpu", train_cfg.get("batch_size", 1)))
        grad_accum = int(ds_cfg.get("gradient_accumulation_steps", train_cfg.get("gradient_accumulation_steps", 1)))
        if micro_batch != int(train_cfg.get("batch_size", micro_batch)):
            raise ValueError("DeepSpeed train_micro_batch_size_per_gpu must match training.batch_size.")
        if grad_accum != int(train_cfg.get("gradient_accumulation_steps", grad_accum)):
            raise ValueError("DeepSpeed gradient_accumulation_steps must match training.gradient_accumulation_steps.")

        ds_config = {
            "train_micro_batch_size_per_gpu": micro_batch,
            "gradient_accumulation_steps": grad_accum,
            "gradient_clipping": float(ds_cfg.get("gradient_clipping", train_cfg.get("gradient_clip_norm", 1.0))),
            "zero_optimization": {
                "stage": int(ds_cfg.get("zero_stage", 2)),
                "offload_optimizer": {"device": "cpu" if ds_cfg.get("offload_optimizer", False) else "none"},
                "offload_param": {"device": "cpu" if ds_cfg.get("offload_param", False) else "none"},
            },
            "fp16": {
                "enabled": bool(ds_cfg.get("fp16_enabled", True)),
                "initial_scale_power": int(ds_cfg.get("initial_scale_power", 16)),
                "loss_scale_window": int(ds_cfg.get("loss_scale_window", 1000)),
                "hysteresis": int(ds_cfg.get("hysteresis", 2)),
                "min_loss_scale": int(ds_cfg.get("min_loss_scale", 1)),
            },
            "bf16": {"enabled": bool(ds_cfg.get("bf16_enabled", False))},
            "torch_autocast": {"enabled": bool(train_cfg.get("amp_enabled", torch.cuda.is_available()))},
        }
        self.model, self.optimizer, _, _ = deepspeed.initialize(
            model=self.model,
            optimizer=self.optimizer,
            config=ds_config,
        )
        self.ds_engine = self.model
        self.scaler = GradScaler("cuda", enabled=False)
        self.logger.info("DeepSpeed ZeRO-%d initialized", ds_config["zero_optimization"]["stage"])

    @staticmethod
    def _schedule_epoch(epoch: int) -> int:
        return max(int(epoch) - 1, 0)

    def _apply_warmup_lr(self, schedule_epoch: int) -> None:
        if self.scheduler_warmup_epochs <= 0 or schedule_epoch >= self.scheduler_warmup_epochs:
            return
        if self.scheduler_warmup_epochs == 1:
            factor = 1.0
        else:
            progress = float(schedule_epoch) / float(self.scheduler_warmup_epochs - 1)
            factor = self.scheduler_warmup_start_factor + (1.0 - self.scheduler_warmup_start_factor) * progress
        for group, base_lr in zip(self.optimizer.param_groups, self._base_lrs):
            group["lr"] = base_lr * factor

    def _step_scheduler(self, val_loss: float) -> None:
        if self.scheduler is None:
            return
        if self._active_scheduler_epoch < self.scheduler_warmup_epochs:
            return
        if isinstance(self.scheduler, ReduceLROnPlateau):
            self.scheduler.step(val_loss)
        else:
            self.scheduler.step()

    def _on_epoch_start(self, epoch: int) -> None:
        schedule_epoch = self._schedule_epoch(epoch)
        self._active_scheduler_epoch = schedule_epoch
        self._apply_warmup_lr(schedule_epoch)
        self.criterion.apply_epoch_schedule(schedule_epoch)
        corruption_factor = get_schedule_factor(self.latent_corruption_cfg, schedule_epoch)
        self.current_corruption_sigma = (
            self.base_corruption_sigma * corruption_factor if self.corruption_enabled else 0.0
        )
        if self.local_rank == 0:
            weights = " ".join(f"{key}={value:.4g}" for key, value in self.criterion.current_weights().items())
            lr = self.optimizer.param_groups[0]["lr"] if self.optimizer.param_groups else 0.0
            self.logger.info(
                "Epoch %03d schedules: lr=%.6g corruption_sigma=%.4g %s",
                epoch,
                lr,
                self.current_corruption_sigma,
                weights,
            )

    def resume_from_checkpoint(self, checkpoint_path: str | Path) -> dict[str, Any]:
        ckpt = torch.load(checkpoint_path, map_location=self.device, weights_only=False)
        state = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt))
        state = {key.removeprefix("module."): value for key, value in state.items()}
        missing, unexpected = unwrap_model(self.model).load_state_dict(state, strict=False)
        if missing or unexpected:
            raise RuntimeError(
                f"Resume checkpoint model-state mismatch: missing={missing[:8]} unexpected={unexpected[:8]}"
            )
        self.best_val = float(ckpt.get("best_val_loss", self.best_val))
        self.bad_epochs = int(ckpt.get("bad_epochs", 0))
        self.history = list(ckpt.get("history", self.history))
        self.global_step = int(ckpt.get("global_step", self.global_step))
        epoch = int(ckpt.get("epoch", 0))
        self.start_epoch = epoch + 1
        if self.ds_engine is not None:
            metadata = ckpt.get("metadata", {})
            ds_dir = ckpt.get("deepspeed_checkpoint_dir")
            ds_tag = ckpt.get("deepspeed_checkpoint_tag")
            if isinstance(metadata, dict):
                ds_dir = ds_dir or metadata.get("deepspeed_checkpoint_dir")
                ds_tag = ds_tag or metadata.get("deepspeed_checkpoint_tag")
            if ds_dir and ds_tag and Path(ds_dir).exists():
                load_path, _ = self.ds_engine.load_checkpoint(
                    load_dir=str(ds_dir),
                    tag=str(ds_tag),
                    load_optimizer_states=True,
                    load_lr_scheduler_states=True,
                    load_module_only=False,
                )
                if self.local_rank == 0:
                    self.logger.info("Resumed DeepSpeed native state from %s tag=%s load_path=%s", ds_dir, ds_tag, load_path)
            elif self.local_rank == 0:
                self.logger.warning("DeepSpeed resume requested but native shard metadata was unavailable; model weights restored only.")
        else:
            if "optimizer_state_dict" in ckpt:
                self.optimizer.load_state_dict(ckpt["optimizer_state_dict"])
            scheduler_state = ckpt.get("scheduler_state_dict")
            if scheduler_state and self.scheduler is not None:
                self.scheduler.load_state_dict(scheduler_state)
        ema_state = ckpt.get("ema_state_dict")
        if self.ema is not None:
            if isinstance(ema_state, dict) and ema_state:
                self.ema.load_state_dict(ema_state, unwrap_model(self.model))
                if self.local_rank == 0:
                    self.logger.info("Restored EMA state from checkpoint.")
            else:
                self.ema.reset_to_model(unwrap_model(self.model))
                if self.local_rank == 0:
                    self.logger.warning("Resume checkpoint has no EMA state; EMA reset to resumed model weights.")
        if self.local_rank == 0:
            self.logger.info(
                "Resumed portable checkpoint %s from epoch %d; missing=%d unexpected=%d",
                checkpoint_path,
                epoch,
                len(missing),
                len(unexpected),
            )
        return ckpt

    def _amp_context(self):
        return autocast(device_type="cuda", dtype=self.amp_dtype, enabled=self.amp_enabled and self.device.type == "cuda")

    def _move_batch(self, batch) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None, torch.Tensor | None]:
        if len(batch) == 4:
            atm, ocn, atm_mask, ocn_mask = batch
        else:
            atm, ocn = batch
            atm_mask = None
            ocn_mask = None
        atm = atm.to(self.device, non_blocking=True)
        ocn = ocn.to(self.device, non_blocking=True)
        if atm_mask is not None:
            atm_mask = atm_mask.to(self.device, non_blocking=True)
        if ocn_mask is not None:
            ocn_mask = ocn_mask.to(self.device, non_blocking=True)
        if getattr(self, "channels_last", False):
            atm = atm.contiguous(memory_format=torch.channels_last)
            ocn = ocn.contiguous(memory_format=torch.channels_last)
            if atm_mask is not None:
                atm_mask = atm_mask.contiguous(memory_format=torch.channels_last)
            if ocn_mask is not None:
                ocn_mask = ocn_mask.contiguous(memory_format=torch.channels_last)
        return atm, ocn, atm_mask, ocn_mask

    def _forward_losses(
        self,
        atm: torch.Tensor,
        ocn: torch.Tensor,
        train: bool,
        atm_mask: torch.Tensor | None = None,
        ocn_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        corruption_sigma = self.current_corruption_sigma if train else 0.0
        single_pass = (
            train
            and corruption_sigma > 0.0
            and self.corruption_mode in {"single_pass", "single_pass_denoising", "denoising"}
        )
        out = self.model(
            atm,
            ocn,
            corruption_sigma=corruption_sigma,
            single_pass_denoising=single_pass,
        )
        return self.criterion(
            out["atm_recon"],
            out["ocn_recon"],
            atm,
            ocn,
            None if single_pass else out.get("atm_corrupt_recon"),
            None if single_pass else out.get("ocn_corrupt_recon"),
            atm_mask=atm_mask,
            ocn_mask=ocn_mask,
        )

    def train_epoch(self, epoch: int) -> dict[str, float]:
        self.model.train()
        if hasattr(self.train_loader.sampler, "set_epoch"):
            self.train_loader.sampler.set_epoch(epoch)
        totals: dict[str, torch.Tensor] = {}
        count = 0
        total_steps = len(self.train_loader)
        epoch_start = time.perf_counter()
        data_total = 0.0
        compute_total = 0.0
        forward_total = 0.0
        backward_total = 0.0
        last_log_step = 0
        last_log_data_total = 0.0
        last_log_compute_total = 0.0
        last_log_forward_total = 0.0
        last_log_backward_total = 0.0
        previous_end = time.perf_counter()
        if self.local_rank == 0:
            lr = self.optimizer.param_groups[0]["lr"] if self.optimizer.param_groups else 0.0
            self.logger.info("Epoch %03d train start: steps=%d lr=%.6g", epoch, total_steps, lr)
        if self.ds_engine is None:
            self.optimizer.zero_grad(set_to_none=True)
        for step, batch in enumerate(self.train_loader, start=1):
            batch_start = time.perf_counter()
            data_wait = batch_start - previous_end
            data_total += data_wait
            if data_wait > self.slow_data_threshold_s:
                self.logger.warning(
                    "Slow DataLoader wait: E%03d [%4d/%4d] rank=%d waited=%.1fs threshold=%.1fs",
                    epoch,
                    step,
                    total_steps,
                    self.local_rank,
                    data_wait,
                    self.slow_data_threshold_s,
                )
            atm, ocn, atm_mask, ocn_mask = self._move_batch(batch)
            compute_start = time.perf_counter()
            with self._amp_context():
                losses = self._forward_losses(atm, ocn, train=True, atm_mask=atm_mask, ocn_mask=ocn_mask)
                loss = losses["total_loss"]
            forward_elapsed = time.perf_counter() - compute_start
            backward_start = time.perf_counter()
            if self.ds_engine is not None:
                self.ds_engine.backward(loss)
                is_update_step = bool(
                    self.ds_engine.is_gradient_accumulation_boundary()
                    if hasattr(self.ds_engine, "is_gradient_accumulation_boundary")
                    else True
                )
                self.ds_engine.step()
                if is_update_step:
                    self.global_step += 1
                    if self.ema is not None:
                        self.ema.update(unwrap_model(self.model))
            else:
                is_update_step = step % self.gradient_accumulation_steps == 0 or step == total_steps
                sync_context = (
                    self.ddp_model.no_sync()
                    if self.ddp_model is not None and not is_update_step
                    else nullcontext()
                )
                with sync_context:
                    self.scaler.scale(loss / self.gradient_accumulation_steps).backward()
                if is_update_step:
                    if self.gradient_clip_norm > 0:
                        self.scaler.unscale_(self.optimizer)
                        torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.gradient_clip_norm)
                    self.scaler.step(self.optimizer)
                    self.scaler.update()
                    self.optimizer.zero_grad(set_to_none=True)
                    self.global_step += 1
                    if self.ema is not None:
                        self.ema.update(unwrap_model(self.model))
            backward_elapsed = time.perf_counter() - backward_start
            forward_total += forward_elapsed
            backward_total += backward_elapsed
            for key, value in losses.items():
                detached = value.detach().float()
                totals[key] = totals.get(key, detached.new_tensor(0.0)) + detached
            count += 1
            compute_elapsed = time.perf_counter() - compute_start
            compute_total += compute_elapsed
            if compute_elapsed > self.slow_compute_threshold_s:
                self.logger.warning(
                    "Slow train compute/sync: E%03d [%4d/%4d] rank=%d elapsed=%.1fs threshold=%.1fs",
                    epoch,
                    step,
                    total_steps,
                    self.local_rank,
                    compute_elapsed,
                    self.slow_compute_threshold_s,
                )
            should_log = (
                step == 1
                or step == total_steps
                or (is_update_step and self.global_step % max(self.log_every_n_steps, 1) == 0)
            )
            if should_log and self.local_rank == 0:
                loss_values = {key: float(value.detach().cpu()) for key, value in losses.items()}
                component_text = " ".join(
                    f"{key}={value:.4g}"
                    for key, value in loss_values.items()
                    if key != "total_loss"
                    and (
                        (key.endswith("_weighted_loss") and not key.startswith("corrupt_"))
                        or key.endswith("_total")
                        or key == "corrupt_total"
                    )
                    and (key == "corrupt_total" or not key.startswith("corrupt_"))
                )
                lr = self.optimizer.param_groups[0]["lr"] if self.optimizer.param_groups else 0.0
                data_ms = data_total / max(step, 1) * 1000.0
                compute_ms = compute_total / max(step, 1) * 1000.0
                window_steps = max(step - last_log_step, 1)
                window_data_ms = (data_total - last_log_data_total) / window_steps * 1000.0
                window_compute_ms = (compute_total - last_log_compute_total) / window_steps * 1000.0
                window_forward_ms = (forward_total - last_log_forward_total) / window_steps * 1000.0
                window_backward_ms = (backward_total - last_log_backward_total) / window_steps * 1000.0
                bottleneck = "DATA" if data_ms > compute_ms else "COMPUTE"
                self.logger.info(
                    "E%03d [%4d/%4d] S%d | total=%.4g %s | lr=%.6g data=%.1fms compute=%.1fms window_data=%.1fms window_compute=%.1fms window_forward=%.1fms window_backward=%.1fms [%s]",
                    epoch,
                    step,
                    total_steps,
                    self.global_step,
                    loss_values.get("total_loss", float(loss.detach().cpu())),
                    component_text,
                    lr,
                    data_ms,
                    compute_ms,
                    window_data_ms,
                    window_compute_ms,
                    window_forward_ms,
                    window_backward_ms,
                    bottleneck,
                )
                last_log_step = step
                last_log_data_total = data_total
                last_log_compute_total = compute_total
                last_log_forward_total = forward_total
                last_log_backward_total = backward_total
            previous_end = time.perf_counter()
        if self.local_rank == 0 and count > 0:
            avg_tensor = totals.get("total_loss")
            avg = float(avg_tensor.cpu()) / max(count, 1) if avg_tensor is not None else 0.0
            self.logger.info(
                "Epoch %03d train done: avg_total=%.6g steps=%d elapsed=%.1fs",
                epoch,
                avg,
                count,
                time.perf_counter() - epoch_start,
            )
        return {key: float(value.cpu()) / max(count, 1) for key, value in totals.items()}

    @torch.no_grad()
    def validate(self) -> dict[str, float]:
        self.model.eval()
        totals: dict[str, torch.Tensor] = {}
        count = 0
        skipped = 0
        start = time.perf_counter()
        if self.local_rank == 0:
            self.logger.info("Validation start: steps=%d", len(self.val_loader))
        for step, batch in enumerate(self.val_loader, start=1):
            atm, ocn, atm_mask, ocn_mask = self._move_batch(batch)
            with self._amp_context():
                losses = self._forward_losses(atm, ocn, train=False, atm_mask=atm_mask, ocn_mask=ocn_mask)
            nonfinite = [key for key, value in losses.items() if not bool(torch.isfinite(value.detach()).all())]
            if nonfinite:
                skipped += 1
                if self.local_rank == 0 and skipped <= 8:
                    self.logger.warning(
                        "Validation non-finite loss skipped: step=%d/%d keys=%s",
                        step,
                        len(self.val_loader),
                        ",".join(nonfinite[:8]),
                    )
                continue
            for key, value in losses.items():
                detached = value.detach().float()
                totals[key] = totals.get(key, detached.new_tensor(0.0)) + detached
            count += 1
        if self.local_rank == 0 and count > 0:
            avg_tensor = totals.get("total_loss")
            avg = float(avg_tensor.cpu()) / max(count, 1) if avg_tensor is not None else 0.0
            self.logger.info(
                "Validation done: avg_total=%.6g steps=%d skipped=%d elapsed=%.1fs",
                avg,
                count,
                skipped,
                time.perf_counter() - start,
            )
        elif self.local_rank == 0:
            self.logger.warning(
                "Validation produced no finite batches: steps=%d skipped=%d elapsed=%.1fs",
                len(self.val_loader),
                skipped,
                time.perf_counter() - start,
            )
        return {key: float(value.cpu()) / max(count, 1) for key, value in totals.items()}

    def _record_and_plot(self, epoch: int, train_metrics: dict[str, float], val_metrics: dict[str, float]) -> None:
        row = {f"train_{k}": v for k, v in train_metrics.items()}
        row.update({f"val_{k}": v for k, v in val_metrics.items()})
        row["total_loss"] = val_metrics.get("total_loss", train_metrics.get("total_loss", float("inf")))
        self.history.append(row)
        self.plotter.plot(self.history, filename=f"losses_epoch{epoch}.png")
        self.plotter.plot(self.history, filename="losses.png")

    def _write_run_files(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        if self.local_rank != 0:
            return
        self.dataset.export_metadata(self.output_dir / "normalization_metadata.json")
        with open(self.output_dir / "config.json", "w", encoding="utf-8") as f:
            json.dump(self.cfg, f, indent=2)

    def _save_checkpoint_with_ema(
        self,
        path: Path,
        epoch: int,
        extra: dict[str, Any],
    ) -> None:
        raw_model = unwrap_model(self.model)
        ema_state = self.ema.state_dict() if self.ema is not None else None
        if ema_state is not None:
            self.ema.apply_shadow(raw_model)
        try:
            payload_extra = dict(extra)
            if self.scheduler is not None:
                payload_extra["scheduler_state_dict"] = self.scheduler.state_dict()
            if ema_state is not None:
                payload_extra["ema_state_dict"] = ema_state
            save_checkpoint(path, self.model, self.optimizer, epoch, self.best_val, extra=payload_extra)
        finally:
            if ema_state is not None:
                self.ema.restore(raw_model)

    @staticmethod
    def _checkpoint_epoch_sort_key(path: Path) -> tuple[int, int, str]:
        stem = path.stem
        prefix = "checkpoint_epoch"
        if stem.startswith(prefix):
            suffix = stem[len(prefix):]
            if suffix.isdigit():
                return (0, int(suffix), path.name)
        return (1, len(path.name), path.name)

    def _cleanup_old_checkpoints(self) -> None:
        if self.local_rank != 0 or self.keep_last_n <= 0:
            return
        ckpts = sorted(
            self.checkpoint_dir.glob("checkpoint_epoch*.pth"),
            key=self._checkpoint_epoch_sort_key,
        )
        if len(ckpts) > self.keep_last_n:
            for old in ckpts[:-self.keep_last_n]:
                old.unlink(missing_ok=True)

        ds_dir = self.checkpoint_dir / "deepspeed"
        if ds_dir.exists():
            ds_ckpts = sorted(
                (path for path in ds_dir.glob("epoch_*") if path.is_dir()),
                key=lambda path: self._checkpoint_epoch_sort_key(Path(path.name.replace("epoch_", "checkpoint_epoch") + ".pth")),
            )
            if len(ds_ckpts) > self.keep_last_n:
                import shutil

                for old in ds_ckpts[:-self.keep_last_n]:
                    shutil.rmtree(old, ignore_errors=True)

    def train(self) -> Path:
        self._write_run_files()
        if self.local_rank == 0:
            self.logger.info("Training loop starting at epoch=%d ending at epoch=%d", self.start_epoch, self.epochs)
        for epoch in range(self.start_epoch, self.epochs + 1):
            self._on_epoch_start(epoch)
            train_metrics = self.train_epoch(epoch)
            self._release_epoch_memory("post-train", epoch)
            should_validate = epoch % self.val_every_n_epochs == 0 or epoch == self.epochs
            if should_validate:
                val_metrics = self.validate()
                self._release_epoch_memory("post-val", epoch)
                val_loss = val_metrics.get("total_loss", float("inf"))
                self._step_scheduler(val_loss)
                is_best = val_loss < self.best_val - self.early_min_delta
                if is_best:
                    self.best_val = val_loss
                    self.bad_epochs = 0
                else:
                    self.bad_epochs += 1
            else:
                val_metrics = {}
                val_loss = None
                is_best = False
                if self.local_rank == 0:
                    self.logger.info(
                        "Validation skipped at epoch %03d: val_every_n_epochs=%d",
                        epoch,
                        self.val_every_n_epochs,
                    )
            deepspeed_meta = {}
            if self.ds_engine is not None:
                ds_dir = self.checkpoint_dir / "deepspeed"
                ds_tag = f"epoch_{epoch:04d}"
                self.ds_engine.save_checkpoint(
                    str(ds_dir),
                    tag=ds_tag,
                    client_state={
                        "epoch": int(epoch),
                        "best_val_loss": float(self.best_val),
                        "bad_epochs": int(self.bad_epochs),
                        "global_step": int(self.global_step),
                        "history": self.history,
                    },
                )
                deepspeed_meta = {
                    "deepspeed_checkpoint_dir": str(ds_dir),
                    "deepspeed_checkpoint_tag": ds_tag,
                }
            if self.local_rank == 0:
                self._record_and_plot(epoch, train_metrics, val_metrics)
                checkpoint_extra = {
                    "history": self.history,
                    "bad_epochs": self.bad_epochs,
                    "global_step": self.global_step,
                    **deepspeed_meta,
                }
                self._save_checkpoint_with_ema(self.checkpoint_dir / "last_model.pth", epoch, checkpoint_extra)
                if epoch % self.save_every_n_epochs == 0:
                    self._save_checkpoint_with_ema(
                        self.checkpoint_dir / f"checkpoint_epoch{epoch}.pth",
                        epoch,
                        checkpoint_extra,
                    )
            best_deepspeed_meta = {}
            if is_best and self.ds_engine is not None:
                ds_dir = self.checkpoint_dir / "deepspeed"
                self.ds_engine.save_checkpoint(
                    str(ds_dir),
                    tag="best",
                    client_state={
                        "epoch": int(epoch),
                        "best_val_loss": float(self.best_val),
                        "bad_epochs": int(self.bad_epochs),
                        "global_step": int(self.global_step),
                        "history": self.history,
                    },
                )
                best_deepspeed_meta = {
                    "deepspeed_checkpoint_dir": str(ds_dir),
                    "deepspeed_checkpoint_tag": "best",
                }
            if is_best and self.local_rank == 0:
                checkpoint_extra = {
                    "history": self.history,
                    "bad_epochs": self.bad_epochs,
                    "global_step": self.global_step,
                    **best_deepspeed_meta,
                }
                self._save_checkpoint_with_ema(self.checkpoint_dir / "best_model.pth", epoch, checkpoint_extra)
            if self.local_rank == 0:
                val_text = "skipped" if val_loss is None else f"{val_loss:.6g}"
                best_text = f"{self.best_val:.6g}" if math.isfinite(self.best_val) else "inf"
                self.logger.info("epoch=%d train=%.6g val=%s best=%s", epoch, train_metrics["total_loss"], val_text, best_text)
                if epoch % self.save_every_n_epochs == 0:
                    self._cleanup_old_checkpoints()
            self._release_epoch_memory("post-checkpoint", epoch)
            if should_validate and self.early_stopping_enabled and self.bad_epochs >= self.early_patience:
                self.logger.info("Early stopping at epoch %d", epoch)
                break
        return self.output_dir
