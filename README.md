# Hf2p_CM2.1_SCDA & WRF_LDA

This repository contains the core code and minimal supporting materials used in the paper *Python–Fortran Hybrid Programming for Deep Incorporation of AI and Physics Modeling and Data Assimilation*. It includes the trimmed and modified files necessary to understand and reproduce parts of the work. **Note:** the observation datasets, pretrained model weights, and the full source trees of CM2.1 and WRF v3.7.1 are **not** included — only the code files we changed are provided here.

To reproduce the experiments you will need to prepare the observation data, train the models locally, and adapt the Python dataflow to your environment.

---

## Associated paper and authors

**Python–Fortran Hybrid Programming for Deep Incorporation of AI and Physics Modeling and Data Assimilation**

by Xianrui Zhu, Zikuan Lin, Shaoqing Zhang, Zebin Lu, Songhua Wu, Xiangyun Hou, Zhisheng Xiao, Zhicheng Ren, Jingyu Li, Jing Xu, Yang Gao, Rixu Hao, Xiaolin Yu, Mingkui Li

---

## Contents

This repository provides minimal implementation files and driver scripts referenced in the paper for:

* `CM2.1_SCDA`
* `WRF_LDA`

Each folder contains the pared-down / modified code and example drivers required to run the manuscript’s demo cases.

Additionally, each of the above folders contains two environment reference files (one Conda and one pip) that can be used to configure a Python environment:

* **For `cm2.1_scda`**

  * `environment_cm2_scda.yaml` (Conda)
  * `requirements.txt` (pip)
* **For `WRF_lda`**

  * `environment_wrflda.yaml` (Conda)
  * `requirements.txt` (pip)

## Important files & layout

* `WRFv3.7.1/cm2.1-modified-src` trimmed/modified CM2.1/WRF code.
* `PMC_w_LDA(SCDA)/` — Python modules and scripts (see notes below).
* `plug/` F2py-related Fortran wrappers used to build `.so` modules (for example `plug.F90` and the corresponding generated shared objects).
* `data/` — partial data used to illustrate/assist the Python implementations (not the full observation sets).

## Notes on Python code and data

* The Python implementations used in the experiments are located in the `PMC_w_LDA(SCDA)/` folder. This includes:

  * Fortran–Python interaction code and the Python main controller.
  * Implementation of the VAE-LDA algorithm used in the manuscript.
  * Other ML model code used in experiments.
* The `data/` folder contains limited example files to help understand the Python-side implementation and dataflow. **These are not the full observation datasets** — full observations must be prepared separately.
* When you run the Python code you will likely need to adapt the file paths and the data flow to match your local filesystem and data preparation pipeline.

## Build / compilation notes

### For `CM2.1_SCDA`

* **Important:** copy the files from `cm2.1-modified-src` (in this repository) into the corresponding locations in your WRF source trees, **overwriting** the original source files, and then compile. The modified files must replace the source files prior to building.
* Build and compile using the **Intel** compiler toolchain.
* Install and compile NetCDF (and the NetCDF–Python interface), OpenMPI, and other compiled dependencies with the Intel compilers to ensure binary compatibility with the provided Fortran code and libraries.
* Avoid relying exclusively on Conda-provided binaries for compiled components — the Intel toolchain is required.
* The compilation/compile-flow is described in Section 2 of the paper — follow that flow for best reproducibility.

### For `WRF_LDA`

* When running WRF’s `./configure` script, **select option `34 1`**. We cannot guarantee correct behavior if a different configuration option is chosen.
* **Important:** copy the files from `WRFv3.7.1-modified-src` (in this repository) into the corresponding locations in your WRF source trees, **overwriting** the original source files, and then compile. The modified files must replace the source files prior to building.
* The compilation/compile-flow is described in Section 2 of the paper — follow that flow for best reproducibility.

## Contact

If you have questions or encounter issues, please contact:

* Shaoqing Zhang — [szhang@ouc.edu.cn](mailto:szhang@ouc.edu.cn)
* Xianrui Zhu — [zhuxianrui@stu.ouc.edu.cn](mailto:zhuxianrui@stu.ouc.edu.cn) / [mapzhu@foxmail.com](mailto:mapzhu@foxmail.com)
