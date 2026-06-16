# CM2-LDA

Lightweight coupled atmosphere-ocean deterministic AE for the copied `temp/cm2_cda_main.py` DA workflow. The project is intentionally small: no score-based DA variants, no VAE path, and no large PALM-GDA code copy.

## Layout

Only these files live in the project root:

- `cm2_cda_main.py`: DA runtime main entry copied from `temp/cm2_cda_main.py`.
- `vae_cda.py`: compatibility adapter exposing `latent_space_da.do_da(...)`.
- `README.md`: this file.

Subpackages:

- `configs/`: path, preprocessing, training, and NMC configs.
- `data_preprocess/`: NetCDF/observation preprocessing and NMC background statistics.
- `model/`: coupled AE, ConvNeXt/residual blocks, cross-attention, and 2D RoPE.
- `training/`: AE trainer, checkpointing, logging, and reconstruction/spectral losses.
- `analysis/`: reconstruction plots, spectra, and metric evaluation.
- `utils/`: path, data, normalization, plotting, and logging utilities.

## Linux Environment

Use the same Python environment style as `palm-GDA`. DeepSpeed is optional for single-GPU runs but required for the DeepSpeed launch below.

```bash
conda activate palm-gda
cd /path/to/VAE/CM2-lda
```

If DeepSpeed needs to build CUDA extensions on the Linux training machine, make sure `CUDA_HOME` points at the CUDA toolkit used by the PALM-GDA environment.

## Path Configuration

All project paths are set in `configs/paths.json`, not in the shell command. Default storage is under `/data/cm2_lda`:

- raw CM2/observation files: `/data/cm2_lda/raw`
- preprocessed frame dataset: `/data/cm2_lda/data`
- default AE checkpoint for analysis/DA: `/data/cm2_lda/checkpoints/run_ae_3x/best_model.pth`
- AE checkpoints: `/data/cm2_lda/checkpoints/run_ae_3x` and `/data/cm2_lda/checkpoints/run_ae_8x`
- NMC output: `/data/cm2_lda/nmc_outputs`
- default analysis output: under the selected checkpoint run directory, for example `/data/cm2_lda/checkpoints/run_ae_3x/analysis_result`

Edit `configs/paths.json`, `configs/preprocess_nc.json`, and `configs/preprocess_observation.json` if your Linux raw-data layout differs.

## Data Preparation

Preprocess CM2 NetCDF into one ordered `.pt` file per frame:

```bash
python -m data_preprocess.preprocess_nc \
  --config configs/preprocess_nc.json
```

Recommended MPI preprocessing for large frame datasets:

```bash
mkdir -p logs
nohup mpirun -np 8 env PYTHONUNBUFFERED=1 python -u -m data_preprocess.preprocess_nc \
  --config configs/preprocess_nc.json \
  --parallel-backend mpi \
  > logs/log.preprocess_nc_mpi8.nohup 2>&1 &
```

MPI launchers may still buffer aggregated child-rank stdout until process exit or kill. Treat `logs/log.preprocess_nc_mpi8.nohup` as launcher diagnostics only. The preprocessor writes reliable per-rank progress files under `logs/`:

```bash
tail -f logs/preprocess_nc_rank00.log
tail -f logs/preprocess_nc_rank*.log
```

If MPI is not available, use local multiprocessing:

```bash
mkdir -p logs
nohup python -m data_preprocess.preprocess_nc \
  --config configs/preprocess_nc.json \
  --parallel-backend process \
  --preprocess-workers 8 \
  > logs/log.preprocess_nc_process8.nohup 2>&1 &
```

The preprocessing follows a two-pass streaming design in parallel mode. Pass 1 lets each MPI rank or local worker read its assigned NetCDF files, compute z-score statistics, and return local file/frame metadata in the same pass. Rank 0 then merges statistics and orders the gathered metadata by source-file index to assign deterministic global frame IDs. Pass 2 rereads chunks, applies the saved z-score, writes `atm_mask`/`ocn_mask`, and saves one normalized `frame_XXXXXXXX.pt` per time step. Frame filenames are assigned by global ordered index, so outputs remain deterministic across ranks. It writes `normalization_stats.json` as a separate sidecar before frame writing starts, while `normalization_metadata.json` and `_frame_manifest.json` only reference that sidecar. No large aggregate `.pt` file is produced. `_frame_manifest.json` records `frame_files`, `time_index`, and `time_values`; NMC uses this order, so do not shuffle or rename frame files.

Preprocess observations for the DA adapter:

```bash
python -m data_preprocess.preprocess_observation \
  --config configs/preprocess_observation.json
```

This writes `/data/cm2_lda/data/merged_observation.pt` with the same `{"data": {"ps1": ..., "sst": ...}}` layout used by `temp/vae_cda.py`.

## AE Training

Use `configs/train_ae_compression_3x.yaml` as the default production training profile. It is the lower-compression config and keeps logs/checkpoints separate from the 8x comparison run.

Recommended long Linux run, 4 GPUs, DeepSpeed launcher with DDP backend, nohup:

```bash
mkdir -p logs
nohup deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_3x.yaml \
  > logs/log.train_ae_3x_4gpu.deepspeed 2>&1 &
```

Follow progress:

```bash
tail -f logs/log.train_ae_3x_4gpu.deepspeed
```

The nohup log is the primary training log. It includes launcher output, config/data/model startup messages, epoch start/end summaries, validation summaries, and training step lines every `logging.log_every_n_steps` steps. The default interval is `50`; each epoch also logs the first and final training step.

Resume the recommended 3x 4-GPU run from the last portable checkpoint:

```bash
mkdir -p logs
nohup deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_3x.yaml \
  --resume-from /data/cm2_lda/checkpoints/run_ae_3x/last_model.pth \
  > logs/log.train_ae_3x_4gpu_resume.deepspeed 2>&1 &
```

The resume mechanism follows the PALM-GDA pattern: `checkpoint.resume_from` points at the canonical rank-0 checkpoint, usually `last_model.pth` or `best_model.pth`. On the DDP path, model, optimizer, scheduler, EMA, epoch, step, and history are restored from the portable rank-0 checkpoint.

Checkpoint saving also follows the PALM-GDA pattern. Each run writes:

- `last_model.pth`: latest epoch, overwritten every epoch.
- `best_model.pth`: best validation epoch so far.
- `checkpoint_epochN.pth`: periodic epoch snapshots controlled by `checkpoint.save_every_n_epochs`; current AE configs save every epoch and keep the latest `checkpoint.keep_last_n` snapshots.

Use an epoch snapshot when comparing spectra at the same training point, for example:

```bash
python -m analysis.analyze_ae \
  --checkpoint /data/cm2_lda/checkpoints/run_ae_3x_cnnsubband_v1/checkpoint_epoch10.pth
```

Compression-profile configs:

- `configs/train_ae_compression_8x.yaml`: `latent_channels=36`, latent dim `40,500`, overall compression about `8.39:1`.
- `configs/train_ae_compression_3x.yaml`: `latent_channels=96`, latent dim `108,000`, overall compression about `3.15:1`.
- `configs/train_ae_compression_3x_cnnsubband.yaml`: code name `cm2_ae_3x_cnnsubband_v1`; same 3x latent size, replaces cross-attention fusion with CNN fusion and uses PALM-style spatial subband loss for high-frequency comparison.
- `configs/train_ae.yaml`: convenience default, currently mirrors the 8x architecture but writes to the generic `run_ae` paths.

The 3x production config defaults to:

- `data.path: /data/cm2_lda/data`
- checkpoint output: `/data/cm2_lda/checkpoints/run_ae_3x`
- decoder-heavy coupled AE: `latent_channels=96`, latent grid `25x45`, coupling grid `50x90`, encoder `width=56`, decoder `width=160`, decoder middle depth `6`, decoder stage depth `[8, 6, 4, 2, 2, 2]`
- `distributed.backend: ddp`
- DeepSpeed is used as the launcher; training uses PyTorch DDP unless `distributed.backend` is explicitly changed to `deepspeed`
- FP16 autocast with DDP gradient accumulation; DDP uses `no_sync()` between accumulation boundaries to avoid unnecessary per-micro-step all-reduce
- `backends.attention: auto`, which tries FlashAttention, then xFormers, then PyTorch SDPA
- `backends.conv: super_conv`, which enables dense Conv2d channels-last optimization for this regular-grid AE

Reference launch commands:

Single GPU, 3x config, no distributed launcher:

```bash
python -m training.train_ae \
  --config configs/train_ae_compression_3x.yaml
```

Because no launcher sets `WORLD_SIZE`, this falls back to normal single-GPU training even though the config says `ddp`.

DDP, 2 GPUs, 3x config through the DeepSpeed launcher:

```bash
deepspeed --num_gpus=2 --module training.train_ae \
  --config configs/train_ae_compression_3x.yaml
```

DDP, 4 GPUs, foreground, 3x config:

```bash
deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_3x.yaml
```

DDP, 4 GPUs, nohup, 8x comparison:

```bash
mkdir -p logs
nohup deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_8x.yaml \
  > logs/log.train_ae_8x_4gpu.deepspeed 2>&1 &
```

DDP, 4 GPUs, nohup, 3x CNN-fusion/spatial-subband comparison:

```bash
mkdir -p logs
nohup deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_3x_cnnsubband.yaml \
  > logs/log.train_ae_3x_cnnsubband_v1_4gpu.deepspeed 2>&1 &
```

Resume the 8x comparison run:

```bash
mkdir -p logs
nohup deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_8x.yaml \
  --resume-from /data/cm2_lda/checkpoints/run_ae_8x/last_model.pth \
  > logs/log.train_ae_8x_4gpu_resume.deepspeed 2>&1 &
```

Equivalent config-file resume, if you prefer editing YAML instead of passing `--resume-from`:

```yaml
checkpoint:
  save_dir: /data/cm2_lda/checkpoints/run_ae_3x
  resume_from: /data/cm2_lda/checkpoints/run_ae_3x/last_model.pth
```

Single GPU, generic default config:

```bash
python -m training.train_ae \
  --config configs/train_ae.yaml
```

DDP, 4 GPUs, generic default config:

```bash
deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae.yaml
```

DeepSpeed launcher is the default multi-GPU entrypoint. With `distributed.backend: ddp`, it launches one process per GPU while training uses PyTorch DDP:

```bash
deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_3x.yaml
```

If you switch a config to `distributed.backend: deepspeed`, DeepSpeed native checkpoint shards are written under the selected checkpoint directory, for example:

```bash
/data/cm2_lda/checkpoints/run_ae_3x/deepspeed/
```

The portable PyTorch checkpoints `best_model.pth` and `last_model.pth` are still saved on rank 0 for analysis, NMC, and DA runtime loading.

## Runtime Backends

The backend block follows the PALM-GDA pattern: it is a set of runtime hints, not mandatory dependencies.

```yaml
runtime:
  allow_tf32: true
  cudnn_benchmark: true
  float32_matmul_precision: high
  channels_last: true
  torch_compile: false
  torch_compile_mode: reduce-overhead
  log_backend_status: true

backends:
  attention: auto
  conv: super_conv
  sparse: spconv
```

Attention backend choices:

- `auto`: try `flash_attn`, then `xformers`, then PyTorch `scaled_dot_product_attention`.
- `flash_attn`: use FlashAttention when installed and compatible with the current dtype/device.
- `xformers`: use `xformers.ops.memory_efficient_attention` when available.
- `sdpa`: use PyTorch SDPA only.

Convolution backend choices:

- `super_conv` or `channels_last`: use channels-last memory format, TF32, and cuDNN benchmark for dense 2D Conv/ConvNeXt blocks.
- `torch`: keep standard PyTorch memory format.
- `runtime.torch_compile: true`: additionally compile the AE with `torch.compile`.

`backends.sparse: spconv` is kept for PALM-GDA config parity. CM2-LDA currently operates on dense 2D atmosphere/ocean grids, so sparse convolution backends are not used unless the model/data representation is changed to sparse tensors.

Training writes:

- `best_model.pth`
- `last_model.pth`
- DeepSpeed native checkpoints under `deepspeed/` only when `distributed.backend: deepspeed`
- `config.json`
- `normalization_metadata.json`
- loss plots under the selected config log directory, for example `logs/run_ae_3x/loss_plots`
- the primary nohup log under local `./logs/`, for example `logs/log.train_ae_3x_4gpu.deepspeed`

## NMC Background

After training, estimate latent-space NMC background statistics:

```bash
python -m data_preprocess.nmc_background \
  --config configs/nmc_background.json
```

The NMC script defaults to `split="val"` with `shuffle=False`, then computes lagged perturbations in strict time order. It saves the full latent `B_diag` only; local covariance blocks are only for figures unless `save_local_plot_block` or `save_diagonal_blocks` is enabled. Off-diagonal diagnostics use the same NMC perturbations but sample many global latent-dimension pairs, controlled by `diagnostic_pair_count`, so changing `plot_dims` changes only the displayed local block and not the global p95 statistic. For full-diagonal visual checks, `diagonal_block_size` and `diagonal_block_count` control tiled covariance/correlation blocks sampled along the full flattened latent diagonal; the default uses 24 evenly spaced 2048x2048 blocks to keep the local view wide while avoiding overly sparse middle coverage.

## Analysis

Per-variable original/reconstruction/difference plots and spectra. By default, outputs are written next to the selected checkpoint under `<run_dir>/analysis_result`:

```bash
python -m analysis.analyze_ae
```

The analysis entrypoint reads `configs/analysis.json` by default. Use that file to select validation-only analysis, frame ranges, and plotting colors. For example, set `"sample_start": 40`, `"sample_stop": 100`, and `"num_samples": null` to analyze validation frames 40 through 99.

Default output layout:

- plots: `/data/cm2_lda/checkpoints/run_ae_3x/analysis_result/plots/reconstruction`
- spectra: `/data/cm2_lda/checkpoints/run_ae_3x/analysis_result/plots/spectra`
- reconstruction metrics: `/data/cm2_lda/checkpoints/run_ae_3x/analysis_result/metrics/reconstruction`
- spectral band metrics: `/data/cm2_lda/checkpoints/run_ae_3x/analysis_result/metrics/spectral_bands`

Use `--output-dir` only when you want to override the run-local analysis folder:

```bash
python -m analysis.analyze_ae \
  --output-dir /data/cm2_lda/analysis_outputs/manual_check
```

Metric-only reconstruction evaluation:

```bash
python -m analysis.evaluate_reconstruction
```

This writes to `<run_dir>/analysis_result/metrics/evaluation` by default.

The evaluator reports MAE, median AE, p95 AE, MSE, RMSE, centered RMSE, normalized errors, bias, SMAPE, MAPE, max absolute error, relative L2, correlation, R2, explained variance, NSE, KGE, Willmott d, global SSIM, PSNR, and spectral L1.

## DA Runtime

Run the copied CM2 DA entry from this directory:

```bash
python cm2_cda_main.py
```

The DA adapter loads paths through `configs/paths.json`: observation `.pt`, AE checkpoint, normalization metadata, and NMC covariance. It keeps the `latent_space_da` interface expected by `cm2_cda_main.py`.
