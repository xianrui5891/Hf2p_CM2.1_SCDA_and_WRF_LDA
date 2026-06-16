import os
import re
import torch
import torch.nn.functional as F
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime, timedelta

from networks.transformer import LGUnet_all
from tqdm import tqdm

# ==== 标准化参数 ====
norm_dict = {
    'U10': {'mean': -0.7289, 'std': 2.0134},
    'V10': {'mean': -0.4459, 'std': 2.2311},
    'T2':  {'mean': 288.7260, 'std': 10.4469},
    'Q2':  {'mean': 0.0101, 'std': 0.0050}
}

# ==== 加载变量 ====
def load_nc_variable(file_path, variable_name='U10'):
    nc_file = Dataset(file_path, 'r')
    if variable_name not in nc_file.variables:
        nc_file.close()
        raise KeyError(f"{variable_name} not found in {file_path}")
    data = nc_file.variables[variable_name][:]
    nc_file.close()

    # 如果是带掩膜的 MaskedArray，先填充（用 NaN 或者 0，根据需求）
    if isinstance(data, np.ma.MaskedArray):
        data = data.filled(np.nan)   # 把掩膜位置设为 NaN（随后用 nan_to_num 处理）

    # 处理时间维度（如果是 3D: time,x,y）
    if data.ndim == 3:
        data = data[0]  # 使用第 1 个时间点

    # 如果 shape 是 (nx, ny) 或 (ny, nx)，保持为 2D ndarray
    # 去掉最后一行/列（你原代码的逻辑）
    data = np.asarray(data)  # 确保是 ndarray
    if data.shape[0] > 0 and data.shape[1] > 0:
        data = data[:-1, :-1]

    return data   # 返回普通的 numpy ndarray，可能包含 NaN

def save_visualization(output_dir_base, filename, var_list, original_inputs, output_denorm, save_fig=True):
    """
    output_dir_base: 基础输出目录，比如 'output_figs'
    filename: 当前.nc文件名，用于命名文件
    var_list: 变量列表，比如 ['U10','V10','T2','Q2']
    original_inputs: list，每个变量的输入图 (H,W)
    output_denorm: list，每个变量的重建图 (H,W)
    """

    for i, var in enumerate(var_list):
        # 创建变量对应的子文件夹
        output_dir = os.path.join(output_dir_base, var)
        os.makedirs(output_dir, exist_ok=True)

        plt.figure(figsize=(10, 4))
        plt.suptitle(f"{filename} - {var}")

        vmin = np.min(original_inputs[i])
        vmax = np.max(original_inputs[i])

        plt.subplot(1, 2, 1)
        plt.imshow(original_inputs[i], cmap='viridis', vmin=vmin, vmax=vmax)
        plt.title(f'Input {var}')
        plt.colorbar()

        plt.subplot(1, 2, 2)
        plt.imshow(output_denorm[i], cmap='viridis', vmin=vmin, vmax=vmax)
        plt.title(f'Reconstructed {var}')
        plt.colorbar()

        # 文件名中去掉.nc后缀，保证文件名干净
        base_name = os.path.splitext(filename)[0]
        save_path = os.path.join(output_dir, f"{base_name}_{var}.png")
        plt.savefig(save_path)
        plt.close()

# ==== Resize 输入 ====
def resize_input(data, target_size=(256, 256)):
    """
    输入可以是：
      - numpy ndarray, shape (H,W) 或 (H,W) 含 NaN
      - torch.Tensor, shape (H,W) 或 (1,H,W) 或 (N,C,H,W)
    返回：
      - torch.Tensor, shape (H,W)    （便于后续统一处理）
    """
    # 如果是 numpy，先转换并保证 dtype float32
    if isinstance(data, np.ndarray):
        # 将 NaN 替换为 0（或使用更合适的填充值），避免插值时抛错
        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
        t = torch.from_numpy(data.astype(np.float32))
    elif isinstance(data, torch.Tensor):
        t = data.to(dtype=torch.float32)
    else:
        raise TypeError("resize_input expects numpy ndarray or torch.Tensor")

    # 现在 t 的 shape 应是 (H,W) 或 (1,H,W) 或 (C,H,W) 或 (N,C,H,W)
    if t.dim() == 2:
        t = t.unsqueeze(0).unsqueeze(0)   # -> (1,1,H,W)
        squeezed_after = True
    elif t.dim() == 3:
        # 可能是 (C,H,W) 或 (1,H,W)；我们期望 (1,C,H,W)
        t = t.unsqueeze(0)               # -> (1,C,H,W)
        squeezed_after = False
    elif t.dim() == 4:
        # 已经是 (N,C,H,W)
        squeezed_after = False
    else:
        raise ValueError(f"Unsupported tensor dim: {t.dim()}")

    # 使用 interpolate 进行 resize（双线性）
    t_resized = F.interpolate(t, size=target_size, mode='bilinear', align_corners=False)

    # 去掉 batch 维度返回 (C,H,W) 或 (1,H,W) 以便后续使用
    t_resized = t_resized.squeeze(0)  # -> (C,H,W) 或 (1,H,W)

    # 如果原始是单通道，返回 shape (H,W) 更方便（但后续代码中需要注意）
    if t_resized.shape[0] == 1:
        return t_resized.squeeze(0)  # -> (H,W)
    else:
        return t_resized            # -> (C,H,W)

# ==== 构建模型 ====
def build_model():
    model = LGUnet_all(
        img_size=[256, 256],
        patch_size=4,
        stride=[4, 4],
        in_chans=4,
        out_chans=4,
        enc_depths=[2, 2, 6],
        enc_heads=[3, 6, 12],
        lg_depths=[2, 2],
        lg_heads=[6, 12],
        inchans_list=[4],
        outchans_list=[4],
        enc_dim=96,
        embed_dim=768,
        window_size=8,
        Weather_T=1,
        use_checkpoint=False,
        pre_norm=True
    )
    return model

# ==== 推理与可视化 ====
def test_model_on_file(model, nc_path, var_list=['U10', 'V10', 'T2', 'Q2'], save_fig=False, output_dir='output_figs'):
    input_channels = []
    original_inputs = []

    for var in var_list:
        raw = load_nc_variable(nc_path, var)   # numpy ndarray
        # 将 NaN 替换（load 已填 NaN），再做 resize（resize 会把 NaN->0）
        # raw = np.nan_to_num(raw)   # 已在 resize 中处理，这里可以不重复
        raw_resized = resize_input(raw, target_size=(256,256))  # 返回 torch.Tensor (H,W)

        # 保存原始（用于可视化）。把 tensor 转回 numpy：注意可能是 torch.Tensor
        original_inputs.append(raw_resized.detach().cpu().numpy() if isinstance(raw_resized, torch.Tensor) else np.asarray(raw_resized))

        # 归一化：先把 tensor 变为 float32 并在通道维上处理
        if isinstance(raw_resized, np.ndarray):
            t = torch.from_numpy(raw_resized.astype(np.float32))
        else:
            t = raw_resized.to(dtype=torch.float32)

        # 现在 t shape 为 (H,W)，把它变成 (C,H,W) 即 (1,H,W)
        t = t.unsqueeze(0)

        mean = norm_dict[var]['mean']
        std = norm_dict[var]['std']
        t = (t - mean) / std   # 广播生效

        input_channels.append(t)   # t shape: (1,H,W)

    # 把所有通道按 channel 维合并 -> (C,H,W)
    input_tensor = torch.cat(input_channels, dim=0)   # (4, H, W)
    # 加 batch 维 -> (1, C, H, W)
    input_tensor = input_tensor.unsqueeze(0).to(dtype=torch.float32)

    # 推理
    model.eval()
    with torch.no_grad():
        output = model(input_tensor).detach().cpu().numpy()[0]  # shape: [C, H, W]

    # 反归一化（convert to numpy）
    output_denorm = []
    for i, var in enumerate(var_list):
        std = norm_dict[var]['std']
        mean = norm_dict[var]['mean']
        recovered = output[i] * std + mean  # numpy运算
        output_denorm.append(recovered)

    # 保存图像，替换原有的显示/保存逻辑
    if save_fig:
        save_visualization(output_dir, os.path.basename(nc_path), var_list, original_inputs, output_denorm, save_fig)
    else:
        for i, var in enumerate(var_list):
            plt.figure(figsize=(10, 4))
            plt.suptitle(f"{os.path.basename(nc_path)} - {var}")

            vmin = np.min(original_inputs[i])
            vmax = np.max(original_inputs[i])

            plt.subplot(1, 2, 1)
            plt.imshow(original_inputs[i], cmap='viridis', vmin=vmin, vmax=vmax)
            plt.title(f'Input {var}')
            plt.colorbar()

            plt.subplot(1, 2, 2)
            plt.imshow(output_denorm[i], cmap='viridis', vmin=vmin, vmax=vmax)
            plt.title(f'Reconstructed {var}')
            plt.colorbar()

            plt.show()
            plt.close()

# ==== 主逻辑 ====
def main():
    model_path = 'model_forgoal.pth'
    data_dir = '/data/zsq_data/d01_d02_d03/d03_wrf'
    start_date = datetime(2024, 9, 30)
    end_date = datetime(2024, 10, 1)
    var_list = ['U10', 'V10', 'T2', 'Q2']

    # 加载模型
    model = build_model()
    checkpoint = torch.load(model_path, map_location='cpu')
    if any(k.startswith('module.') for k in checkpoint.keys()):
        from collections import OrderedDict
        new_state_dict = OrderedDict()
        for k, v in checkpoint.items():
            new_state_dict[k[7:]] = v
        checkpoint = new_state_dict
    model.load_state_dict(checkpoint)

    # 遍历文件夹下所有 .nc 文件
    all_files = sorted(os.listdir(data_dir))
    pattern = re.compile(r'(\d{4}-\d{2}-\d{2})')

    for fname in tqdm(all_files, desc='Processing files'):
        match = pattern.search(fname)
        if not match:
            continue

        file_date = datetime.strptime(match.group(1), "%Y-%m-%d")
        if start_date <= file_date <= end_date:
            full_path = os.path.join(data_dir, fname)
            try:
                test_model_on_file(model, full_path, var_list=var_list, save_fig=True)
            except Exception as e:
                print(f"❌ Failed on {fname}: {e}")

if __name__ == '__main__':
    main()