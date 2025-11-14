# Hf2p_CM2.1_SCDA & WRF_LDA

This repository contains the core code and minimal supporting materials used in the paper *Python–Fortran Hybrid Programming for Deep Incorporation of AI and Physics Modeling and Data Assimilation*. It includes the trimmed and modified files necessary to understand and reproduce parts of the work. **Note:** the observation datasets, pretrained model weights, and the full source trees of CM2.1 and WRF v3.7.1 are **not** included — only the code files we changed are provided here.

To reproduce the experiments you will need to prepare the observation data, train the models locally, and adapt the Python dataflow to your environment.

## Contents

* Minimal implementation files and driver scripts referenced in the paper for:

  * `CM2.1_SCDA`
  * `WRF_LDA`

Each folder contains the pared-down/modified code and example drivers required to run the manuscript’s demo cases.

## Prerequisites / Important notes

### For `CM2.1_SCDA`

* Build and compile with the **Intel** compiler toolchain.
* Install and compile NetCDF (and the NetCDF–Python interface), OpenMPI, and other dependencies using the Intel compilers.
* Avoid relying exclusively on Conda-provided binaries for compiled components — the Intel toolchain is required to ensure binary compatibility with the provided code.

### For `WRF_LDA`

* When running WRF’s `./configure` script, select option **`34 1`**. We cannot guarantee correct behavior if a different configuration option is chosen.

## Reproducibility

This repository provides only the modified code and minimal scripts used in the study. To fully reproduce the experiments you must:

1. Acquire or prepare the observation datasets used in the paper.
2. Install and build the required software using the toolchain and options described above.
3. Train the AI models locally (pretrained weights are not provided).

## Contact

If you have questions or encounter issues, please contact:

* Shaoqing Zhang — [szhang@ouc.edu.cn](mailto:szhang@ouc.edu.cn)
* Xianrui Zhu — [zhuxianrui@stu.ouc.edu.cn](mailto:zhuxianrui@stu.ouc.edu.cn) / [mapzhu@foxmail.com](mailto:mapzhu@foxmail.com)

---

Thank you for your interest — contributions, bug reports, and questions are welcome.
