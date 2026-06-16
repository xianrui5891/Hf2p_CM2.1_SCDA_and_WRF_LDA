import os
import numpy as np
import netCDF4 as nc
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset
from datetime import datetime


class WRFSlidingDataset(Dataset):
    def __init__(self, nc_dir, input_len=1, label_len=1,
                 input_vars=None, label_vars=None,
                 start_month=1, end_month=12):

        self.nc_dir = nc_dir
        self.input_len = input_len
        self.label_len = label_len
        self.input_vars = input_vars
        self.label_vars = label_vars
        self.target_height = 256
        self.target_width = 256

        self.nc_files = []
        self.sample_start_times = []  # 每个样本的起始时间

        # === 收集符合时间范围的文件 ===
        all_files = sorted(os.listdir(nc_dir))
        for fname in all_files:
            try:
                time_str = fname.split("_")[-1]
                date_str = fname.split("_")[-2]
                dt = datetime.strptime(date_str + "_" + time_str, "%Y-%m-%d_%H:%M:%S")
                if start_month <= dt.month <= end_month:
                    self.nc_files.append(os.path.join(nc_dir, fname))
            except Exception as e:
                print(f"[跳过] 无法解析时间: {fname}, 错误: {e}")

        # === 计算样本数 ===
        self.total_samples = len(self.nc_files) - self.input_len - self.label_len + 1
        assert self.total_samples > 0, f"季度筛选后文件不足，仅有 {len(self.nc_files)} 个"

        # === 记录每个样本的起始时间 ===
        for i in range(self.total_samples):
            fname = os.path.basename(self.nc_files[i])
            date_str = fname.split("_")[-2]
            time_str = fname.split("_")[-1]
            dt = datetime.strptime(date_str + "_" + time_str, "%Y-%m-%d_%H:%M:%S")
            self.sample_start_times.append(dt)

        print(f"[INFO] 数据加载完成，季度内样本数: {self.total_samples}")

    def __len__(self):
        return self.total_samples

    def __getitem__(self, idx):
        input_seq = []

        for i in range(self.input_len):
            input_seq.append(self._load_vars(self.nc_files[idx + i], self.input_vars))

        # concat -> shape [C_in * input_len, H, W]
        input_tensor = torch.from_numpy(np.concatenate(input_seq, axis=0)).float()

        # === 使用双线性插值调整到 256x256（替代原来的 padding 逻辑） ===
        # interpolate 需要输入为 4D: (N, C, H, W)
        input_tensor = input_tensor.unsqueeze(0)  # (1, C, H, W)
        input_tensor = F.interpolate(input_tensor, size=(self.target_height, self.target_width), mode='bilinear', align_corners=False)
        input_tensor = input_tensor.squeeze(0)  # (C, H, W)

        label_tensor = input_tensor.clone()

        # 使用插值后没有填充区域，mask 全为 1
        mask = torch.ones_like(label_tensor)

        return input_tensor.float(), label_tensor.float(), mask.float()

    def _load_vars(self, filepath, vars_list):
        stats = {
            'U10': {'mean': -0.7289, 'std': 2.0134},
            'V10': {'mean': -0.4459, 'std': 2.2311},
            'T2': {'mean': 288.7260, 'std': 10.4469},
            'Q2': {'mean': 0.0101, 'std': 0.0050}
        }

        with nc.Dataset(filepath) as ds:
            var_data = []
            for var in vars_list:
                data = ds[var][:]
                # 如果是 3 维（例如 time, y, x），取第一个 time 层
                if data.ndim == 3:
                    # 只取 data[0] 作为二维场
                    data = data[0, :, :]
                elif data.ndim == 2:
                    # 要求：取全部值 -1 即不取最后一行和最后一列
                    data = data[:-1, :-1]
                else:
                    raise ValueError(f"[维度异常] {filepath} 中变量 {var} 维度为 {data.shape}")

                data = np.nan_to_num(data)

                mean = stats[var]['mean']
                std = stats[var]['std']
                data = (data - mean) / std

                var_data.append(data)
        return np.stack(var_data, axis=0)

    def get_sample_quarters(self, base_year=2022):
        """
        返回每个样本所属的季度编号（从 base_year 开始，每3个月为一个季度）
        """
        quarters = []
        for dt in self.sample_start_times:
            delta = (dt.year - base_year) * 12 + (dt.month - 1)
            quarter_id = delta // 3 + 1  # 每三个月为一个季度
            quarters.append(quarter_id)
        return quarters
