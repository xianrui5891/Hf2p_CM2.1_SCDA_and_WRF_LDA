# Hf2p CM2.1-SCDA and WRF-LDA

This repository contains the hybrid Python-Fortran code used for CM2.1-SCDA and WRF-LDA experiments in the Hf2pMDA work. The CM2.1 Python DA component has been replaced with the CM2-LDA autoencoder workflow from `D:\Desktop\models\VAE\CM2-lda`, including preprocessing, training, analysis, NMC background estimation, runtime DA adapter, and configs.

Observation data, processed training frames, normalization metadata, trained model checkpoints, and NMC `B_z`/`.npz` files are prepared by the provided scripts and are not stored in this repository. The public data record is available at Zenodo: https://doi.org/10.5281/zenodo.18799861.

## Paper

**Python-Fortran Hybrid Programming for Deep Incorporation of AI and Physics Modeling and Data Assimilation (Hf2pMDA_1.0)**

Authors: Xianrui Zhu, Zikuan Lin, Shaoqing Zhang, Zebin Lu, Songhua Wu, Xiangyun Hou, Zhisheng Xiao, Zhicheng Ren, Jiangyu Li, Jing Xu, Yang Gao, Rixu Hao, Xiaolin Yu, Mingkui Li.

## Repository Layout

- `CM2.1-SCDA/`: CM2.1-SCDA source modifications, F2PY plug files, Python CM2-LDA workflow, and CM2 environment references.
- `CM2.1-SCDA/cm2.1-modified-src/`: modified CM2.1 Fortran source tree used before building CM2.1-SCDA.
- `CM2.1-SCDA/plug/`: F2PY interface files for the CM2.1 Python-Fortran bridge.
- `CM2.1-SCDA/PMC_w_SCDA/`: current CM2-LDA Python workflow and runtime DA adapter.
- `WRF_LDA/`: WRF v3.7.1 source modifications, F2PY plug files, Python LDA workflow, and WRF environment references.
- `WRF_LDA/WRFv3.7.1-modified-src/`: modified WRF source files.
- `WRF_LDA/PMC_w_LDA/`: WRF-LDA Python workflow.

Environment references:

- `CM2.1-SCDA/environment_cm2_scda.yaml`
- `CM2.1-SCDA/requirements.txt`
- `WRF_LDA/environment_wrflda.yaml`
- `WRF_LDA/requirements.txt`

## CM2-LDA Python Workflow

The current CM2.1 Python workflow lives in `CM2.1-SCDA/PMC_w_SCDA`. It contains:

- `cm2_cda_main.py`: main CM2.1-SCDA runtime controller called from the coupled Python-Fortran workflow.
- `vae_cda.py`: DA adapter exposing `latent_space_da.do_da(...)`.
- `configs/`: path, preprocessing, training, analysis, and NMC configuration files.
- `data_preprocess/`: CM2 NetCDF preprocessing, observation preprocessing, NMC background, and climatological background scripts.
- `training/`: autoencoder training entrypoint, trainer, losses, optimizer, and checkpoint utilities.
- `model/`: coupled atmosphere-ocean autoencoder model.
- `analysis/`: reconstruction analysis, metrics, and plotting utilities.
- `utils/`: path, data, logging, normalization, EMA, and plotting helpers.
Generated logs, analysis result folders, plot folders, model checkpoints, normalization metadata, NMC output files, and `__pycache__` files are intentionally not kept here.

## CM2-LDA Paths

Most CM2-LDA scripts read paths from JSON/YAML files under `CM2.1-SCDA/PMC_w_SCDA/configs`. These configs still contain absolute paths from the original Linux training/runtime environment, such as `/data/cm2_lda` and `/data/ouc/...`.

Before running on a new machine, update the path fields in these files:

- `configs/paths.json`
- `configs/preprocess_nc.json`
- `configs/preprocess_observation.json`
- `configs/train_ae_compression_8x.yaml`
- `configs/nmc_background.json`
- `configs/analysis.json`

For runtime DA, `configs/paths.json` controls the observation `.pt`, AE checkpoint, normalization metadata, normalization statistics, and NMC covariance `.npz` paths used by `vae_cda.py`.

Key path fields to check:

- `configs/paths.json`: `storage_root`, `raw_dir`, `data_dir`, `runs_dir`, `analysis_dir`, `nmc_output_dir`, `processed_data_path`, `metadata_path`, `normalization_stats_path`, `observation_path`, `ae_model_path`, `nmc_background_covariance_path`.
- `configs/preprocess_nc.json`: `atm_glob`, `ocn_glob`, `output_dir`, `progress_log_dir`.
- `configs/preprocess_observation.json`: `atm_glob`, `ocn_glob`, `output_dir`.
- `configs/train_ae_compression_8x.yaml`: `data.path`, `data.metadata_path`, `checkpoint.save_dir`, `checkpoint.resume_from`, `logging.output_dir`, `logging.log_dir`, `output_dir`.
- `configs/nmc_background.json`: `data_path`, `metadata_path`, `checkpoint`, `output_dir`.
- `configs/analysis.json`: `data_path`, `metadata_path`, `checkpoint`, `output_dir`.

## CM2-LDA Data Preprocessing

Run commands from `CM2.1-SCDA/PMC_w_SCDA`:

```bash
cd CM2.1-SCDA/PMC_w_SCDA
```

Preprocess CM2 NetCDF output into normalized frame files:

```bash
python -m data_preprocess.preprocess_nc --config configs/preprocess_nc.json
```

For MPI preprocessing:

```bash
mkdir -p logs
mpirun -np 8 python -u -m data_preprocess.preprocess_nc \
  --config configs/preprocess_nc.json \
  --parallel-backend mpi
```

For local multiprocessing:

```bash
python -m data_preprocess.preprocess_nc \
  --config configs/preprocess_nc.json \
  --parallel-backend process \
  --preprocess-workers 8
```

Preprocess observations for the DA adapter:

```bash
python -m data_preprocess.preprocess_observation \
  --config configs/preprocess_observation.json
```

This produces the merged observation file expected by `vae_cda.py`.

## CM2-LDA Training

Train the CM2-LDA model with the 8x compression config. This is the target CM2-LDA model version for this workflow:

```bash
python -m training.train_ae \
  --config configs/train_ae_compression_8x.yaml
```

Multi-GPU launch through DeepSpeed using PyTorch DDP internally:

```bash
deepspeed --num_gpus=4 --module training.train_ae \
  --config configs/train_ae_compression_8x.yaml
```

Resume from a checkpoint:

```bash
python -m training.train_ae \
  --config configs/train_ae_compression_8x.yaml \
  --resume-from /path/to/last_model.pth
```

Other training configs are kept only for comparison or development:

- `configs/train_ae_compression_3x.yaml`
- `configs/train_ae_compression_3x_cnnsubband.yaml`
- `configs/train_ae.yaml`

## CM2-LDA NMC Background

Estimate latent-space NMC background statistics after training:

```bash
python -m data_preprocess.nmc_background \
  --config configs/nmc_background.json
```

The script writes `latent_nmc_background_covariance.npz` and diagnostic figures to the configured output directory.

Climatological background variants are also provided:

```bash
python -m data_preprocess.climatological_background \
  --config configs/climatological_background_paper.json

python -m data_preprocess.state_climatological_background \
  --config configs/state_climatological_background_paper.json
```

## CM2-LDA Analysis

Run reconstruction and spectral analysis:

```bash
python -m analysis.analyze_ae --config configs/analysis.json
```

Metric-only reconstruction evaluation:

```bash
python -m analysis.evaluate_reconstruction
```

Additional plotting helpers are available:

- `plot_cm2_training_output_fields.py`
- `plot_obs_model_error.py`

## CM2-LDA Runtime

The coupled CM2.1 Python runtime entry is:

```bash
python cm2_cda_main.py
```

`cm2_cda_main.py` imports the F2PY-generated `cm2` module and calls `cm2_cda_plugs`, so the CM2.1 Fortran side and F2PY bridge must be built and importable before this command can run. Rank 0 constructs `latent_space_da` from `vae_cda.py`; the adapter loads observation data, normalization metadata, the AE model checkpoint, and the NMC covariance path configured in `configs/paths.json`.

## CM2.1 Build Notes

- Copy or merge files from `CM2.1-SCDA/cm2.1-modified-src` into the corresponding CM2.1 source locations before compiling.
- Build the Fortran model and F2PY plug files with a consistent compiler/MPI/NetCDF stack.
- The original workflow used Intel compilers and Intel-compatible NetCDF/OpenMPI builds.
- Avoid mixing incompatible Conda binary libraries with the Fortran/MPI stack.

## WRF-LDA Notes

- WRF modified source files are under `WRF_LDA/WRFv3.7.1-modified-src`.
- Python workflow files are under `WRF_LDA/PMC_w_LDA`.
- F2PY bridge files are under `WRF_LDA/plug`.
- When configuring WRF v3.7.1, use option `34 1` for the workflow documented in this project.
- Copy or merge files from `WRFv3.7.1-modified-src` into the corresponding WRF source locations before building.

## Data

Large raw and preprocessed observation/model datasets are not tracked here. Use the Zenodo record and then update the config paths for your local filesystem:

https://doi.org/10.5281/zenodo.18799861

## Contact

- Shaoqing Zhang: [szhang@ouc.edu.cn](mailto:szhang@ouc.edu.cn)
- Xianrui Zhu: [zhuxianrui@stu.ouc.edu.cn](mailto:zhuxianrui@stu.ouc.edu.cn), [mapzhu@foxmail.com](mailto:mapzhu@foxmail.com)

## Citation

```bibtex
@Article{egusphere-2025-6479,
AUTHOR = {Zhu, X. and Lin, Z. and Zhang, S. and Lu, Z. and Wu, S. and Hou, X. and Xiao, Z. and Ren, Z. and Li, J. and Xu, J. and Gao, Y. and Hao, R. and Yu, X. and Li, M.},
TITLE = {Python-Fortran Hybrid Programming for Deep Incorporation of AI and Physics Modeling and Data Assimilation (Hf2pMDA\_1.0)},
JOURNAL = {EGUsphere},
VOLUME = {2026},
YEAR = {2026},
PAGES = {1--32},
URL = {https://egusphere.copernicus.org/preprints/2026/egusphere-2025-6479/},
DOI = {10.5194/egusphere-2025-6479}
}
```
