#!/usr/bin/env python3
# 最简脚本：从 .pt 加载 obs_data，逐个 nc_idx 并排画 sst 和 ps1，保存到文件夹

import os
import torch
import numpy as np
import matplotlib.pyplot as plt

obs_path = "/data/ouc/zhuxianrui/CM2-f2py/out_nc_observation/merged_observation.pt"
out_dir = "out_images"
os.makedirs(out_dir, exist_ok=True)

obs = torch.load(obs_path, weights_only=False)
obs_data = obs["data"]

sst_list = obs_data["sst"]
ps1_list = obs_data["ps1"]

def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().squeeze().numpy()
    arr = np.array(x)
    return np.squeeze(arr)

n = len(sst_list)
for idx in range(n):
    sst = to_numpy(sst_list[idx])
    ps1 = to_numpy(ps1_list[idx])

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    im0 = axes[0].imshow(sst, origin="lower")
    axes[0].set_title(f"nc_{idx} sst")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(ps1, origin="lower")
    axes[1].set_title(f"nc_{idx} ps1")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    plt.tight_layout()
    outpath = os.path.join(out_dir, f"nc_{idx:04d}_sst_ps1.png")
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
