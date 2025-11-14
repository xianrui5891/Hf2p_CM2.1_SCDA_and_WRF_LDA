"""
data_output_utils.py

Minimal utilities to save 2D fields into NetCDF files (one file per category: original/da/increment)
and to plot 2D fields.

Modifications requested by user:
- NetCDF: ignore time (no time axis) and do not create per-variable files/folders. Instead create three files
  in the base directory: original.nc, da.nc, increment.nc. Each file contains lon, lat and variables saved by var_name
  (shape: lon x lat). Repeated saves overwrite the variable inside that category file.
- Plotting: use colormap 'viridis'. Make lon the horizontal axis (x) and lat the vertical axis (y).
  Overlayed 'x' markers: black 'x' without annotations.

Assumptions (kept minimal):
- lon_2d[:,0] gives lon axis (1D, length nx)
- lat_2d[0,:] gives lat axis (1D, length ny)
- For 3D fields (u,v,t) shape is (nx, nz, ny) and we plot/save level 0 via u[:,0,:]

Dependencies: numpy, netCDF4, matplotlib

Author: ChatGPT (modified per user request)
"""

import os
import numpy as np
import matplotlib.pyplot as plt


class FieldPlotter:
    """Plot 2D fields with lon as x-axis and lat as y-axis, using colormap 'viridis'.

    Files are saved to: base_dir/<data_type>/<var_name>/<time>_<data_type>_<var_name>.png
    """

    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        os.makedirs(self.base_dir, exist_ok=True)

    def _make_path(self, data_type: str, var_name: str):
        folder = os.path.join(self.base_dir, data_type, var_name)
        os.makedirs(folder, exist_ok=True)
        return folder

    def plot_and_save(self, time_str: str, var_name: str, data2d: np.ndarray,
                    data_type: str = 'original',
                    vmin=None, vmax=None):
        """
        Plot a 2D field and save PNG.

        - lon is horizontal axis (x)
        - lat is vertical axis (y)
        - points: list of (lon_val, lat_val); plotted as black 'x' without annotations.
        - plot_bbox: optional tuple (x1, y1, x2, y2) in lon/lat coordinates to crop the plotted area.
                    Order doesn't matter (function uses min/max).
        Returns saved file path.
        """

        folder = self._make_path(data_type, var_name)
        filename = f"{time_str}_{data_type}_{var_name}.png"
        path = os.path.join(folder, filename)

        arr = np.asarray(data2d)
        if arr.ndim != 2:
            raise ValueError("data2d must be 2D")
        nx, ny = arr.shape

        # 如果需要裁切区域，根据经纬找到索引

        # 转置以满足 imshow(rows=y, cols=x) 的要求
        data_to_plot = arr.T  # shape -> (len(lat_subset), len(lon_subset))

        # extent 的顺序：(x0, x1, y0, y1) —— 可为升序或降序

        fig, ax = plt.subplots(figsize=(7, 5))
        if data_type == 'increment':
            im = ax.imshow(data_to_plot, origin='lower', aspect='auto',
                        vmin=vmin, vmax=vmax, cmap='RdBu_r')
        else: 
            im = ax.imshow(data_to_plot, origin='lower', aspect='auto',
                        vmin=vmin, vmax=vmax, cmap='viridis')
        ax.set_xlabel('lon')
        ax.set_ylabel('lat')
        ax.set_title(f"{time_str} {data_type} {var_name}")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(var_name)

        plt.tight_layout()
        # 确保保存目录存在
        os.makedirs(folder, exist_ok=True)
        fig.savefig(path, dpi=150)
        plt.close(fig)

        return path
