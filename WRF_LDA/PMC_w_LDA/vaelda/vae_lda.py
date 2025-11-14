import torch
import numpy as np
from vaelda.networks.transformer import LGUnet_all
from torch.utils.data import DataLoader
import xarray as xr
from torch.serialization import add_safe_globals
import json
from typing import Union, Dict, Any, Tuple, Optional
add_safe_globals([torch.utils.data.dataset.ConcatDataset])
from utils.normalization_utils import *
from utils.obsdata_utils import *
from utils.interpolate_utils import interpolate as myinterpolate
from utils.debug_utils import DebugPlotter

class latent_space_da:
    def __init__(self, dev = "cuda",
                 learning_rate = 1e-3,
                 max_iterations = 300,
                 obs_path = "./obs/",
                 vae_model_path = "best_model.pth",
                 name_mapping=None):
        
        self.dev = torch.device(dev)
        self.learning_rate = learning_rate
        self.max_iterations = max_iterations
        self.obs_data = obs_data_loader(obs_dir=obs_path).load_all()
        self.var_name = var_name
        self.norm_dict = norm_dict
        self.name_mapping = name_mapping or {'U10':0,'V10':1}
        self.debug_plotter = DebugPlotter(filename_prefix='training_plot')

        self.ckpt = torch.load(vae_model_path, weights_only=False)
        state = self.ckpt
        new_state = {}
        for k, v in state.items():
            if k.startswith("module."):
                new_state[k[len("module."):]] = v
            else:
                new_state[k] = v
        
        self.vae = LGUnet_all(
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
        self.vae.load_state_dict(new_state)
        self.vae.to(self.dev)
        self.vae.eval()
        print("Model loaded.")

    def encode(self, input: torch.tensor):
        input = torch.nn.functional.interpolate(input, size=(256, 256), mode='bilinear', align_corners=False)
        return self.vae.encode(input)

    def decode(self, input: torch.tensor, axis: list):
        decoded = self.vae.decode(input)
        return torch.nn.functional.interpolate(decoded, size=(axis[0].size,axis[1].size), mode='bilinear', align_corners=False)

    @torch.enable_grad()
    def H(self,latent_x, obs_coords, axis: list):
        """
        观测算子 H(latent_x, obs_coords)
        - latent_x: 潜变量（torch.Tensor）
        - obs_coords: [(lat_atm_pts* lon_atm_pts)]
        每个子元素为 1D torch.Tensor，长度为对应观测点数
        返回：
        - 在观测点上的模型值，1D tensor（atm_points）
        """
        
        phy_space = self.decode(latent_x, axis).squeeze(0)  # (C,H,W)
        # 假设 name_mapping 的 ind 是 channel idx
        temp_Hx = []
        for key, ind in self.name_mapping.items():
            channel_field = phy_space[ind]  # (H,W)
            vals = myinterpolate(channel_field, axis, obs_coords, type='bilinear')
            if len(vals)>0:
                temp_Hx.extend(vals)
        Hx = torch.stack(temp_Hx).to(self.dev)  # 1D tensor, requires_grad True if latent_x requires_grad
        print(f"Hx= {Hx.detach().cpu()}")
        print("DEBUG Hx.requires_grad:", Hx.requires_grad)
        return Hx.squeeze()

    @torch.enable_grad()
    def apply_3DVar(self, H, B_diag, R_diag, xb, y, obs_coords, axis, max_iterations=1000, learning_rate=1e-3):
        """
        最小化 J(x) = (x-xb)^T B^{-1} (x-xb) + (y - H(x, obs_coords))^T R^{-1} (y - H(x, obs_coords))
        
        参数：
        H: callable(latent_x, obs_mask) -> obs_vector
        B_diag: 背景协方差矩阵的对角线元素 (state_dim,) - 对角矩阵优化
        R_diag: 观测误差协方差矩阵的对角线元素 (obs_dim,) - 对角矩阵优化
        xb: 背景态（潜空间向量）
        y: 观测值向量（obs_dim,）
        obs_coords: obs值所对应的下标。
        返回：
        最优潜变量（detach clone）
        """
        # 将背景态 xb 复制为可优化参数 x
        x = torch.nn.Parameter(xb.detach().clone())
        optimizer = torch.optim.Adam([x], lr=learning_rate)
        best_x = x.detach().clone()
        B_diag = B_diag.to(x.device)
        R_diag = R_diag.to(x.device)

        print(f"obs: {y.detach().cpu()}")

        for it in range(max_iterations):
            optimizer.zero_grad(set_to_none=True)
            
            # 背景项： (x - xb)^T B^{-1} (x - xb)
            # 对于对角矩阵 B，B^{-1} @ v = v / B_diag
            delta_x = x.ravel() - xb.ravel()
            loss_background = torch.dot(delta_x / B_diag, delta_x)
            
            # 观测项： (y - H(x, obs_mask))^T R^{-1} (y - H(x, obs_mask))
            # 对于对角矩阵 R，R^{-1} @ v = v / R_diag
            model_observations = H(x, obs_coords, axis)
            residual = y.ravel() - model_observations.ravel()
            loss_observation = torch.dot(residual / R_diag, residual)
            
            # 总能量（目标函数）
            loss_J = loss_background + loss_observation
            
            # 反向传播与参数更新
            loss_J.backward(retain_graph=True)
            optimizer.step()
            
            # 保存当前解的副本
            best_x = x.detach().clone()
            
            # 可选：打印调试信息
            if it % 10 == 0:
                print(f"iter {it+1}/{max_iterations}, loss={loss_J.item():.6e}, loss_background={loss_background.item():.6e}, loss_observation = {loss_observation.item():.6e}")
            
            self.debug_plotter.plot(it+1, loss_J.item(), loss_background.item(), loss_observation.item())
        
        return best_x
    
    def to_tensor_from_list(self, input:list) -> torch.Tensor:
        normalized_list = []
        for ind, key in enumerate(self.var_name):
            input[ind] = normalize_standardize(input[ind], self.norm_dict[key]['mean'], self.norm_dict[key]['std'])
            normalized_list.append(input[ind])#注意一下这块！！！可能最后用转置
        
        arr = np.stack(normalized_list, axis=0)
        t = torch.from_numpy(arr).to(torch.float32)
        t = t.unsqueeze(0)  # (1, C, H, W)
        t = t.to(self.dev)
        return t

    def tensor_to_numpy_list(self, input:torch.Tensor) -> list:
        t = input.detach().cpu()
        if t.dim() == 4 and t.shape[0] == 1:
            t = t.squeeze(0)
        denormalized_list = []
        for ind, key in enumerate(self.var_name):
            denormalized = denormalize_standardize(t[ind], self.norm_dict[key]['mean'], self.norm_dict[key]['std'])
            if isinstance(denormalized, torch.Tensor):
                denormalized = denormalized.numpy()

            denormalized_list.append(denormalized.astype(np.float64))#注意一下这块！！！可能最后用转置
        
        return denormalized_list

    @torch.no_grad()
    def do_da(self, time_str:str, input: list, axis: list): 
        """
        对时间帧 nc 做同化
        - nc: 时间索引（会按 self.interval 缩放）
        - input_atm / input_ocn: list of 2D arrays（每个通道一个 (H,W)）
        - axis 格点坐标轴，先x后y，里面是个np.array
        返回：atm_decoded_list, ocn_decoded_list
        """
        # 这种情况不用做同化的，直接退出
        if time_str not in self.obs_data:
            print(f"{time_str} is not in obs_data time set.")
            print(f"obs_data time set:{list(self.obs_data.keys())}")
            return input
        
        # 将输入打包成模型期望的 tensor 形式 (1, C, H, W)
        input_t = self.to_tensor_from_list(input)
        print(f"the input_shape is: {input_t.shape}")
        print(f"the axis[0] shape is: {axis[0].size}")
        print(f"the axis[1] shape is: {axis[1].size}")
        assert input_t.shape[2] == axis[0].size and input_t.shape[3] == axis[1].size

        #记录原本为nan的格点，传回去时要保持原样
        mask_input = ~torch.isnan(input_t).to(self.dev)

        #将nan值赋值为0，和训练时保持一致
        input_t = torch.nan_to_num(input_t, nan=0.0).to(self.dev)

        # 编码获得背景态 xb（潜空间）
        latent_x = self.encode(input_t)

        # 解析前空间大小并获得B_diag
        B_diag = torch.full((latent_x[0].numel(),), 1, device=self.dev)
        
        # 构造观测值和观测值的坐标
        temp_value = {}
        obs_coords = []
        for coor, dicts in self.obs_data[time_str].items():
            y, x = coor
            if axis[0][0] < x < axis[0][-1] and axis[1][0] < y < axis[1][-1]:
                obs_coords.append( (x, y) )
                for key, val in dicts.items():
                    if key not in temp_value:
                        temp_value[key] = []
                    temp_value[key].append(val)
        
        obs_values_list = []
        for key, v in temp_value.items():
            obs_values_list.extend(v)
        obs_values_numpy = np.array(obs_values_list, dtype=np.float32) 
        obs_values = torch.from_numpy(obs_values_numpy).to(self.dev)

        # 构造观测误差协方差矩阵 R（对角小值，保持原始思路）
        R_diag = torch.full((obs_values.shape[0],), 1e-3, device=self.dev)

        # 调用 apply_3DVar
        assim_latent_x = self.apply_3DVar(self.H, B_diag, R_diag, latent_x, obs_values, obs_coords, axis,
                                    self.max_iterations, self.learning_rate)

        self.debug_plotter.close()

        # 解码同化后潜变量S
        decoded = self.decode(assim_latent_x, axis)

        # 重新进行mask
        decoded.masked_fill_(~mask_input, float('nan'))

        # 转回 numpy list 并返回（保持原返回类型）
        decoded_vallist = self.tensor_to_numpy_list(decoded)
        print(f"da had done. the obs_data time set is {time_str}")
        return decoded_vallist
            