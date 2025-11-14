import json
import numpy as np
import torch
from typing import Union, Dict, Any, Tuple, Optional


class NormalizationConfig:
    """归一化配置类，用于加载和管理metadata"""
    
    def __init__(self, metadata_path: Optional[str] = None, metadata_dict: Optional[Dict] = None):
        """
        初始化配置
        
        Args:
            metadata_path: metadata JSON文件路径
            metadata_dict: 或直接传入metadata字典
        """
        if metadata_path:
            with open(metadata_path, 'r') as f:
                self.metadata = json.load(f)
        elif metadata_dict:
            self.metadata = metadata_dict
        else:
            raise ValueError("必须提供 metadata_path 或 metadata_dict")
        
        self.global_min = self.metadata['global_min']
        self.global_max = self.metadata['global_max']
        self.atm_vars = self.metadata['atm_vars']
        self.ocn_vars = self.metadata['ocn_vars']
        self.handle_nan = self.metadata.get('handle_nan', True)
        self.nan_fill_value = self.metadata.get('nan_fill_value', 0.0)
        
    def get_var_range(self, var_name: str, domain: str = 'atm') -> Tuple[float, float]:
        """
        获取变量的最小最大值
        
        Args:
            var_name: 变量名
            domain: 'atm' 或 'ocn'
            
        Returns:
            (min_val, max_val)
        """
        key = f"{domain}::{var_name}"
        if key not in self.global_min or key not in self.global_max:
            raise ValueError(f"变量 {key} 不在metadata中")
        return self.global_min[key], self.global_max[key]


def normalize_minmax(
    data: Union[np.ndarray, torch.Tensor],
    min_val: float,
    max_val: float,
    target_range: Tuple[float, float] = (0.0, 1.0),
    handle_nan: bool = True,
    nan_fill_value: float = 0.0
) -> Union[np.ndarray, torch.Tensor]:
    """
    Min-Max归一化到指定范围
    
    Args:
        data: 输入数据 (H, W)
        min_val: 原始数据最小值
        max_val: 原始数据最大值
        target_range: 目标范围，默认(0, 1)
        handle_nan: 是否处理NaN值
        nan_fill_value: NaN填充值
        
    Returns:
        归一化后的数据
    """
    is_torch = isinstance(data, torch.Tensor)
    
    # 处理NaN
    if handle_nan:
        if is_torch:
            nan_mask = torch.isnan(data)
            data_clean = data.clone()
            data_clean[nan_mask] = nan_fill_value
        else:
            nan_mask = np.isnan(data)
            data_clean = data.copy()
            data_clean[nan_mask] = nan_fill_value
    else:
        data_clean = data
    
    # 归一化
    if max_val - min_val == 0:
        # 避免除零
        normalized = data_clean * 0 + target_range[0]
    else:
        # 标准min-max归一化
        normalized = (data_clean - min_val) / (max_val - min_val)
        # 缩放到目标范围
        target_min, target_max = target_range
        normalized = normalized * (target_max - target_min) + target_min
    
    return normalized


def denormalize_minmax(
    normalized_data: Union[np.ndarray, torch.Tensor],
    min_val: float,
    max_val: float,
    source_range: Tuple[float, float] = (0.0, 1.0)
) -> Union[np.ndarray, torch.Tensor]:
    """
    Min-Max反归一化
    
    Args:
        normalized_data: 归一化的数据 (H, W)
        min_val: 原始数据最小值
        max_val: 原始数据最大值
        source_range: 归一化时的范围，默认(0, 1)
        
    Returns:
        反归一化后的数据
    """
    source_min, source_max = source_range
    
    # 先还原到[0, 1]范围
    if source_max - source_min == 0:
        data_01 = normalized_data * 0
    else:
        data_01 = (normalized_data - source_min) / (source_max - source_min)
    
    # 还原到原始范围
    denormalized = data_01 * (max_val - min_val) + min_val
    
    return denormalized


def normalize_standardize(
    data: Union[np.ndarray, torch.Tensor],
    mean: float,
    std: float,
    handle_nan: bool = True,
    nan_fill_value: float = 0.0
) -> Union[np.ndarray, torch.Tensor]:
    """
    Z-score标准化 (mean=0, std=1)
    
    Args:
        data: 输入数据 (H, W)
        mean: 均值
        std: 标准差
        handle_nan: 是否处理NaN值
        nan_fill_value: NaN填充值
        
    Returns:
        标准化后的数据
    """
    is_torch = isinstance(data, torch.Tensor)
    
    # 处理NaN
    if handle_nan:
        if is_torch:
            nan_mask = torch.isnan(data)
            data_clean = data.clone()
            data_clean[nan_mask] = nan_fill_value
        else:
            nan_mask = np.isnan(data)
            data_clean = data.copy()
            data_clean[nan_mask] = nan_fill_value
    else:
        data_clean = data
    
    # 标准化
    if std == 0:
        standardized = data_clean * 0
    else:
        standardized = (data_clean - mean) / std
    
    return standardized


def denormalize_standardize(
    standardized_data: Union[np.ndarray, torch.Tensor],
    mean: float,
    std: float
) -> Union[np.ndarray, torch.Tensor]:
    """
    Z-score反标准化
    
    Args:
        standardized_data: 标准化的数据 (H, W)
        mean: 均值
        std: 标准差
        
    Returns:
        反标准化后的数据
    """
    denormalized = standardized_data * std + mean
    return denormalized


class VariableNormalizer:
    """基于metadata的变量归一化器"""
    
    def __init__(self, config: NormalizationConfig):
        """
        初始化归一化器
        
        Args:
            config: NormalizationConfig实例
        """
        self.config = config
        
    def normalize(
        self,
        data: Union[np.ndarray, torch.Tensor],
        var_name: str,
        domain: str = 'atm',
        method: str = 'minmax',
        target_range: Tuple[float, float] = (0.0, 1.0)
    ) -> Union[np.ndarray, torch.Tensor]:
        """
        归一化变量
        
        Args:
            data: 输入数据 (H, W)
            var_name: 变量名 (如 'T_SURF', 'SST')
            domain: 'atm' 或 'ocn'
            method: 'minmax' 或 'standardize'
            target_range: 目标范围 (仅用于minmax)
            
        Returns:
            归一化后的数据
        """
        min_val, max_val = self.config.get_var_range(var_name, domain)
        
        if method == 'minmax':
            return normalize_minmax(
                data, min_val, max_val, 
                target_range=target_range,
                handle_nan=self.config.handle_nan,
                nan_fill_value=self.config.nan_fill_value
            )
        elif method == 'standardize':
            # 使用min/max估算mean和std
            mean = (min_val + max_val) / 2
            std = (max_val - min_val) / 6  # 假设6-sigma范围
            return normalize_standardize(
                data, mean, std,
                handle_nan=self.config.handle_nan,
                nan_fill_value=self.config.nan_fill_value
            )
        else:
            raise ValueError(f"不支持的归一化方法: {method}")
    
    def denormalize(
        self,
        normalized_data: Union[np.ndarray, torch.Tensor],
        var_name: str,
        domain: str = 'atm',
        method: str = 'minmax',
        source_range: Tuple[float, float] = (0.0, 1.0)
    ) -> Union[np.ndarray, torch.Tensor]:
        """
        反归一化变量
        
        Args:
            normalized_data: 归一化的数据 (H, W)
            var_name: 变量名
            domain: 'atm' 或 'ocn'
            method: 'minmax' 或 'standardize'
            source_range: 原归一化范围 (仅用于minmax)
            
        Returns:
            反归一化后的数据
        """
        min_val, max_val = self.config.get_var_range(var_name, domain)
        
        if method == 'minmax':
            return denormalize_minmax(
                normalized_data, min_val, max_val,
                source_range=source_range
            )
        elif method == 'standardize':
            mean = (min_val + max_val) / 2
            std = (max_val - min_val) / 6
            return denormalize_standardize(normalized_data, mean, std)
        else:
            raise ValueError(f"不支持的反归一化方法: {method}")


# 使用示例
if __name__ == "__main__":
    # 示例metadata
    metadata = {
        "global_min": {
            "atm::T_SURF": 180.35459899902344,
            "ocn::SST": -1.9007948637008667,
        },
        "global_max": {
            "atm::T_SURF": 351.682373046875,
            "ocn::SST": 35.94911575317383,
        },
        "atm_vars": ["T_SURF"],
        "ocn_vars": ["SST"],
        "handle_nan": True,
        "nan_fill_value": 0.0
    }
    
    # 创建配置和归一化器
    config = NormalizationConfig(metadata_dict=metadata)
    normalizer = VariableNormalizer(config)
    
    # NumPy示例
    data_np = np.random.uniform(200, 300, (90, 144))
    normalized_np = normalizer.normalize(data_np, 'T_SURF', 'atm')
    denormalized_np = normalizer.denormalize(normalized_np, 'T_SURF', 'atm')
    print(f"NumPy - 原始范围: [{data_np.min():.2f}, {data_np.max():.2f}]")
    print(f"NumPy - 归一化范围: [{normalized_np.min():.2f}, {normalized_np.max():.2f}]")
    print(f"NumPy - 反归一化范围: [{denormalized_np.min():.2f}, {denormalized_np.max():.2f}]")
    
    # PyTorch示例
    data_torch = torch.randn(200, 360) * 10 + 15  # SST数据
    normalized_torch = normalizer.normalize(data_torch, 'SST', 'ocn')
    denormalized_torch = normalizer.denormalize(normalized_torch, 'SST', 'ocn')
    print(f"\nTorch - 原始范围: [{data_torch.min():.2f}, {data_torch.max():.2f}]")
    print(f"Torch - 归一化范围: [{normalized_torch.min():.2f}, {normalized_torch.max():.2f}]")
    print(f"Torch - 反归一化范围: [{denormalized_torch.min():.2f}, {denormalized_torch.max():.2f}]")