import torch
import numpy as np
from couple_VAE_large import VAE
from torch.utils.data import DataLoader
from normalization_utils import NormalizationConfig, VariableNormalizer
import xarray as xr
from torch.serialization import add_safe_globals
import json
from typing import Union, Dict, Any, Tuple, Optional
add_safe_globals([torch.utils.data.dataset.ConcatDataset])
from debug_utils import DebugPlotter

class latent_space_da:
    def __init__(self, da_interval : int, dev = "cuda",
                 learning_rate = 1e-3,
                 max_iterations = 30,
                 atm_lat_size = 90,
                 atm_lon_size = 144,
                 ocn_lat_size = 200,
                 ocn_lon_size = 360,
                 obs_path = "../data/merged_observation.pt",
                 vae_model_path = "../data/best_model.pth",
                 metadata_path = "../data/assim_data_metadata.json"):
        
        self.interval = da_interval
        self.dev = torch.device(dev)
        self.learning_rate = learning_rate
        self.max_iterations = max_iterations
        self.obs = torch.load(obs_path, weights_only=False)
        self.obs_data = self.obs["data"]
        # self.obs_coords = self.obs["coords"]
        self.debug_plotter = DebugPlotter(filename_prefix='training_plot')

        with open(metadata_path, 'r', encoding='utf-8') as f:
            self.metadata = json.load(f)
        self.normalizedconfig = NormalizationConfig(metadata_dict=self.metadata)
        self.normalizer = VariableNormalizer(self.normalizedconfig)

        self.atm_vars=self.metadata["atm_vars"]
        self.ocn_vars=self.metadata["ocn_vars"]
    
        # self.lat_atm = self.obs_coords["lat_atm"].to(dev)
        # self.lon_atm = self.obs_coords["lon_atm"].to(dev)
        # self.lat_ocn = self.obs_coords["lat_ocn"].to(dev)
        # self.lon_ocn = self.obs_coords["lon_ocn"].to(dev)
        self.obs_data["ps1"] = self.obs_data["ps1"].squeeze().to(self.dev)
        self.obs_data["sst"] = self.obs_data["sst"].squeeze().to(self.dev)
        self.obs_data["ps1"] = self.normalizer.normalize(self.obs_data["ps1"],
                var_name="PS",
                domain="atm",
                method='minmax')
        self.obs_data["sst"] = self.normalizer.normalize(self.obs_data["sst"],
                var_name="SST",
                domain="ocn",
                method='minmax')
        
        self.ckpt = torch.load(vae_model_path, weights_only=False)
        self.latent_dim = self.ckpt.get("latent_dim", 4096)
        state = self.ckpt.get("model_state_dict", self.ckpt.get("state_dict", None))
        if state is None:
            raise RuntimeError("checkpoint did not contain 'model_state_dict' or 'state_dict'.")
        # strip possible module. prefix
        new_state = {}
        for k, v in state.items():
            if k.startswith("module."):
                new_state[k[len("module."):]] = v
            else:
                new_state[k] = v
        
        self.vae = VAE(latent_dim=self.latent_dim)
        self.vae.load_state_dict(new_state)
        self.vae.to(self.dev)
        self.vae.eval()
        print("Model loaded. latent_dim=", self.latent_dim)

        def H(latent_x, obs_mask):
            """
            观测算子 H(latent_x, obs_coords)
            - latent_x: 潜变量（torch.Tensor）
            - obs_coords: [(lat_atm_pts* lon_atm_pts), (lat_ocn_pts* lon_ocn_pts)]
            每个子元素为 1D torch.Tensor，长度为对应观测点数
            返回：
            - 在观测点上的模型值，1D tensor（atm_points 然后 ocn_points）
            """
           
            obs_atm, obs_ocn = self.vae.decode(latent_x) #(atm/ocn,1, C, H, W)
            
            #print(obs_atm.shape,obs_ocn.shape)
            # 截断到 [0,1]
            obs_atm = obs_atm.clamp(0.0, 1.0).squeeze()#(C, H, W)
            obs_ocn = obs_ocn.clamp(0.0, 1.0).squeeze()#(C, H, W)

            # 构造拟真实值（按相同顺序拼接：atm then ocn）注意，此处数字应跟随metadata位置修改，cm2_cda处也一样
            dev = obs_atm.device
            mask_atm = obs_mask[0].to(dev)
            mask_ocn = obs_mask[1].to(dev)
            Hx_atm = obs_atm[3][mask_atm] #3 为 PS
            Hx_ocn = obs_ocn[0][mask_ocn] #0 为 SST

            # 拼接观测向量（保持 atm 在前，ocn 在后）
            Hx = torch.cat([Hx_atm, Hx_ocn], dim=0)
            return Hx.squeeze()
        
        self.H = H
        #self.B = torch.eye(self.latent_dim, device=dev)
        self.B_diag = torch.full((self.latent_dim,), 1, device=self.dev)
    
    @torch.no_grad()
    def to_tensor_from_list(
        self,
        type: str,
        np_list,
        device: str = "cuda",
        target_range: Tuple[float, float] = (0.0, 1.0)
    ) -> torch.Tensor:
        """
        将numpy数组列表转换为tensor，并进行归一化
        
        Args:
            np_list: numpy数组列表，每个元素形状为(H, W)
            type: 'atm' 或 'ocn'
            normalizer: VariableNormalizer实例
            atm_vars: 大气变量列表
            ocn_vars: 海洋变量列表
            device: 目标设备
            target_range: 归一化目标范围
            
        Returns:
            归一化后的tensor, shape: (1, C, H, W)
        """
        if type not in ['atm', 'ocn']:
            raise ValueError(f"type必须是'atm'或'ocn', 得到: {type}")
        
        # 选择对应的变量列表
        var_list = self.atm_vars if type == 'atm' else self.ocn_vars
        
        if len(np_list) != len(var_list):
            raise ValueError(f"{type}数据通道数({len(np_list)})与变量数({len(var_list)})不匹配")
        
        # 对每个通道进行归一化
        normalized_list = []
        for i, (np_arr, var_name) in enumerate(zip(np_list, var_list)):
            # 归一化当前通道
            normalized = self.normalizer.normalize(
                np_arr,
                var_name=var_name,
                domain=type,
                method='minmax',
                target_range=target_range
            )
            # 转换为numpy（如果是tensor的话）
            if isinstance(normalized, torch.Tensor):
                normalized = normalized.cpu().numpy()
            normalized_list.append(normalized.T)
        
        # Stack成 (C, H, W)
        arr = np.stack(normalized_list, axis=0)
        
        # 转为tensor
        t = torch.from_numpy(arr).to(torch.float32)
        t = t.unsqueeze(0)  # (1, C, H, W)
        
        if device is not None:
            t = t.to(device)
        
        return t

    @torch.no_grad()
    def tensor_to_numpy_list(
        self,
        type: str,
        tensor: torch.Tensor,
        source_range: Tuple[float, float] = (0.0, 1.0)
    ) -> list:
        """
        将tensor转换为numpy数组列表，并进行反归一化
        
        Args:
            type: 'atm' 或 'ocn'
            tensor: 输入tensor, shape: (1, C, H, W) 或 (C, H, W)
            normalizer: VariableNormalizer实例
            atm_vars: 大气变量列表
            ocn_vars: 海洋变量列表
            source_range: 归一化时使用的范围
            
        Returns:
            反归一化后的numpy数组列表，每个元素shape: (H, W), dtype: float64
        """
        if type not in ['atm', 'ocn']:
            raise ValueError(f"type必须是'atm'或'ocn', 得到: {type}")
        
        # 选择对应的变量列表
        var_list = self.atm_vars if type == 'atm' else self.ocn_vars
        
        # 处理tensor维度
        t = tensor.detach().cpu()
        if t.dim() == 4 and t.shape[0] == 1:
            t = t.squeeze(0)  # (C, H, W)
        
        if t.shape[0] != len(var_list):
            raise ValueError(f"{type}数据通道数({t.shape[0]})与变量数({len(var_list)})不匹配")
        
        # 对每个通道进行反归一化
        denormalized_list = []
        for i, var_name in enumerate(var_list):
            # 获取当前通道的tensor
            channel_tensor = t[i]  # (H, W)
            
            # 反归一化
            denormalized = self.normalizer.denormalize(
                channel_tensor,
                var_name=var_name,
                domain=type,
                method='minmax',
                source_range=source_range
            )
            
            # 转换为numpy
            if isinstance(denormalized, torch.Tensor):
                denormalized = denormalized.numpy()
            
            # 转换为float64
            denormalized_list.append(denormalized.astype(np.float64).T)
        
        return denormalized_list
    
    @torch.enable_grad()
    def apply_3DVar(self, H, B_diag, R_diag, xb, y, obs_mask, max_iterations=1000, learning_rate=1e-3):
        """
        最小化 J(x) = (x-xb)^T B^{-1} (x-xb) + (y - H(x, obs_mask))^T R^{-1} (y - H(x, obs_mask))
        
        参数：
        H: callable(latent_x, obs_mask) -> obs_vector
        B_diag: 背景协方差矩阵的对角线元素 (state_dim,) - 对角矩阵优化
        R_diag: 观测误差协方差矩阵的对角线元素 (obs_dim,) - 对角矩阵优化
        xb: 背景态（潜空间向量）
        y: 观测值向量（obs_dim,）
        obs_mask: 经典mask，传给H，构造和y一样大小的，隐空间映射在物理空间的拟真实值的tensor
        返回：
        最优潜变量（detach clone）
        """
        # 将背景态 xb 复制为可优化参数 x
        x = torch.nn.Parameter(xb.detach().clone())
        optimizer = torch.optim.Adam([x], lr=learning_rate)
        best_x = x.detach().clone()
        B_diag = B_diag.to(x.device)
        R_diag = R_diag.to(x.device)

        for it in range(max_iterations):
            optimizer.zero_grad(set_to_none=True)
            
            # 背景项： (x - xb)^T B^{-1} (x - xb)
            # 对于对角矩阵 B，B^{-1} @ v = v / B_diag
            delta_x = x.ravel() - xb.ravel()
            loss_background = torch.dot(delta_x / B_diag, delta_x)
            
            # 观测项： (y - H(x, obs_mask))^T R^{-1} (y - H(x, obs_mask))
            # 对于对角矩阵 R，R^{-1} @ v = v / R_diag
            model_observations = H(x, obs_mask)
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
            # if it % 10 == 0:
            #     print(f"iter {it+1}/{max_iterations}, loss={loss_J.item():.6e}")
            self.debug_plotter.plot(it+1, loss_J.item(), loss_background.item(), loss_observation.item())
        
        return best_x

    @torch.no_grad()
    def do_da(self, nc: int, input_atm: list, input_ocn: list):
        """
        对时间帧 nc 做同化（接口保持不变）
        - nc: 时间索引（会按 self.interval 缩放）
        - input_atm / input_ocn: list of 2D arrays（每个通道一个 (H,W)）
        返回：atm_decoded_list, ocn_decoded_list（与原接口一致）
        """
        # 将输入打包成模型期望的 tensor 形式 (1, C, H, W)
        input_atm_t = self.to_tensor_from_list("atm", input_atm)
        input_ocn_t = self.to_tensor_from_list("ocn", input_ocn)

        #记录原本为nan的格点，传回去时要保持原样
        mask_atm_input = ~torch.isnan(input_atm_t).to(self.dev)
        mask_ocn_input = ~torch.isnan(input_ocn_t).to(self.dev)

        #将nan值赋值为0，和训练时保持一致
        input_atm_t = torch.nan_to_num(input_atm_t, nan=0.0)
        input_ocn_t = torch.nan_to_num(input_ocn_t, nan=0.0)

        # 编码获得背景态 xb（潜空间）
        latent_x, mean, logvar = self.vae.encode(input_atm_t, input_ocn_t)

        # 按 interval 缩放时间索引
        nc_idx = int(nc / self.interval)

        # 生成有效观测点掩码（True 表示有效）
        mask_atm = ~torch.isnan(self.obs_data['ps1'][nc_idx]).to(self.dev)
        mask_ocn = ~torch.isnan(self.obs_data['sst'][nc_idx]).to(self.dev)
        obs_mask = [mask_atm, mask_ocn]

        # # 构造观测坐标（按掩码抽取经纬）
        # lat_grid_atm, lon_grid_atm = torch.meshgrid(self.lat_atm, self.lon_atm, indexing='ij')  # both shape (lat, lon)
        # # 取出有观测的格点对应的 lat/lon（结果是一维向量，长度 = 有效观测数）
        # atm_lat_pts = lat_grid_atm[mask_atm]
        # atm_lon_pts = lon_grid_atm[mask_atm]

        # lat_grid_ocn, lon_grid_ocn = torch.meshgrid(self.lat_ocn, self.lon_ocn, indexing='ij')  # (lat_ocn, lon_ocn)
        # ocn_lat_pts = lat_grid_ocn[mask_ocn]
        # ocn_lon_pts = lon_grid_ocn[mask_ocn]

        # obs_coords = [(atm_lat_pts.to(self.dev), atm_lon_pts.to(self.dev)),
        #         (ocn_lat_pts.to(self.dev), ocn_lon_pts.to(self.dev))]

        # 构造观测值（按相同顺序拼接：atm then ocn）
        obs_atm_vals = self.obs_data['ps1'][nc_idx][mask_atm]
        obs_ocn_vals = self.obs_data['sst'][nc_idx][mask_ocn]
        obs_values = torch.cat([obs_atm_vals, obs_ocn_vals], dim=0).to(self.dev)

        # 构造观测误差协方差矩阵 R（对角小值，保持原始思路）
        R_diag = torch.full((obs_values.shape[0],), 1e-8, device=self.dev)

        # 调用 apply_3DVar
        assim_latent_x = self.apply_3DVar(self.H, self.B_diag, R_diag, latent_x, obs_values, obs_mask,
                                    self.max_iterations, self.learning_rate)

        self.debug_plotter.close()
        # 解码同化后潜变量S
        decoded = self.vae.decode(assim_latent_x)
        atm_decoded_tensor, ocn_decoded_tensor = decoded[0], decoded[1]

        # 重新进行mask
        atm_decoded_tensor.masked_fill_(~mask_atm_input, float('nan'))
        ocn_decoded_tensor.masked_fill_(~mask_ocn_input, float('nan'))

        # 转回 numpy list 并返回（保持原返回类型）
        atm_decoded = self.tensor_to_numpy_list("atm", atm_decoded_tensor)
        ocn_decoded = self.tensor_to_numpy_list("ocn", ocn_decoded_tensor)
        return atm_decoded, ocn_decoded
            