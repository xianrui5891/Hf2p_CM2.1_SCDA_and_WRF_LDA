import numpy as np
from typing import Union, Dict, Any, Tuple, Optional
import torch

var_name = ['U10', 'V10', 'T2', 'Q2']
norm_dict = {
    'U10': {'mean': -0.7289, 'std': 2.0134},
    'V10': {'mean': -0.4459, 'std': 2.2311},
    'T2':  {'mean': 288.7260, 'std': 10.4469},
    'Q2':  {'mean': 0.0101, 'std': 0.0050}
}

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
