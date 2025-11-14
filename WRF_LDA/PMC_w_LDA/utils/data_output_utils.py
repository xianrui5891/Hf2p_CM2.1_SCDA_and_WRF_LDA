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
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


class NCWriter:
    """Write 2D fields into one of three category files (original.nc, da.nc, increment.nc).

    Each category file is at: base_dir/<category>.nc
    Inside each file we store dimensions 'lon' and 'lat' and variables named by var_name with shape (lon, lat).

    Repeated calls with the same var_name will overwrite the variable inside that category file.
    """

    VALID_CATEGORIES = ('original', 'da', 'increment')

    def __init__(self, base_dir: str, lon_2d: np.ndarray, lat_2d: np.ndarray):
        self.base_dir = base_dir
        os.makedirs(self.base_dir, exist_ok=True)

        # extract 1D axes
        self.lon = np.asarray(lon_2d)[:, 0].copy()
        self.lat = np.asarray(lat_2d)[0, :].copy()

        assert self.lon.ndim == 1 and self.lat.ndim == 1, "lon/lat extraction failed"

    def _file_path_for_category(self, category: str):
        assert category in self.VALID_CATEGORIES, f"category must be one of {self.VALID_CATEGORIES}"
        return os.path.join(self.base_dir, f"{category}.nc")

    def save_2d(self, time_str: str, var_name: str, data2d: np.ndarray, data_type: str = "original", units: str = "", long_name: str = ""):
        """
        Save a 2D field into the category file. time_str is ignored for filename (kept for interface compatibility).

        Parameters kept the same for interface compatibility.
        - data_type must be one of: 'original', 'da', 'increment'
        - data2d: shape (nx, ny) where nx == len(self.lon) and ny == len(self.lat)
        """
        category = data_type
        assert category in self.VALID_CATEGORIES, f"data_type must be one of {self.VALID_CATEGORIES}"

        arr = np.asarray(data2d)
        assert arr.ndim == 2, "data2d must be 2D"
        nx, ny = arr.shape
        assert nx == len(self.lon) and ny == len(self.lat), f"data shape {arr.shape} doesn't match lon/lat sizes {(len(self.lon), len(self.lat))}"

        path = self._file_path_for_category(category)

        # If file exists, open in r+ and reuse dimensions/vars; otherwise create new file
        if os.path.exists(path):
            ds = Dataset(path, 'r+')
            try:
                # ensure dims exist
                if 'lon' not in ds.dimensions:
                    ds.createDimension('lon', nx)
                if 'lat' not in ds.dimensions:
                    ds.createDimension('lat', ny)

                # ensure lon/lat vars exist
                if 'lon' not in ds.variables:
                    lonv = ds.createVariable('lon', 'f8', ('lon',))
                    lonv[:] = self.lon
                    lonv.units = 'degrees_east'
                if 'lat' not in ds.variables:
                    latv = ds.createVariable('lat', 'f8', ('lat',))
                    latv[:] = self.lat
                    latv.units = 'degrees_north'

                # create or overwrite variable named var_name with dims (lon, lat)
                if var_name in ds.variables:
                    var = ds.variables[var_name]
                    var[:] = arr
                else:
                    var = ds.createVariable(var_name, 'f4', ('lon', 'lat'))
                    var[:, :] = arr
                    var.units = units
                    var.long_name = long_name

                ds.history = f"Modified {datetime.utcnow().isoformat()}Z"
            finally:
                ds.close()
        else:
            ds = Dataset(path, 'w', format='NETCDF4')
            try:
                ds.createDimension('lon', nx)
                ds.createDimension('lat', ny)

                lonv = ds.createVariable('lon', 'f8', ('lon',))
                latv = ds.createVariable('lat', 'f8', ('lat',))
                lonv[:] = self.lon
                latv[:] = self.lat
                lonv.units = 'degrees_east'
                latv.units = 'degrees_north'

                var = ds.createVariable(var_name, 'f4', ('lon', 'lat'))
                var[:, :] = arr
                var.units = units
                var.long_name = long_name

                ds.history = f"Created {datetime.utcnow().isoformat()}Z"
            finally:
                ds.close()

        return path


class FieldPlotter:
    """Plot 2D fields with lon as x-axis and lat as y-axis, using colormap 'viridis'.

    Files are saved to: base_dir/<data_type>/<var_name>/<time>_<data_type>_<var_name>.png
    """

    def __init__(self, base_dir: str, lon_2d: np.ndarray, lat_2d: np.ndarray):
        self.base_dir = base_dir
        os.makedirs(self.base_dir, exist_ok=True)
        self.lon = np.asarray(lon_2d)[:, 0].copy()
        self.lat = np.asarray(lat_2d)[0, :].copy()

    def _make_path(self, data_type: str, var_name: str):
        folder = os.path.join(self.base_dir, data_type, var_name)
        os.makedirs(folder, exist_ok=True)
        return folder

    def _coord_to_frac_index(self, lon_val, lat_val):
        i_lon = np.interp(lon_val, self.lon, np.arange(len(self.lon)))
        i_lat = np.interp(lat_val, self.lat, np.arange(len(self.lat)))
        return float(i_lon), float(i_lat)

    def plot_and_save(self, time_str: str, var_name: str, data2d: np.ndarray,
                    data_type: str = 'original', points=None,
                    vmin=None, vmax=None, plot_bbox: tuple = None):
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
        if nx != len(self.lon) or ny != len(self.lat):
            raise ValueError("data shape doesn't match lon/lat: "
                            f"data.shape={arr.shape}, len(lon)={len(self.lon)}, len(lat)={len(self.lat)}")

        # 默认使用全域 lon/lat
        lon_arr = np.asarray(self.lon)
        lat_arr = np.asarray(self.lat)

        # 如果需要裁切区域，根据经纬找到索引
        if plot_bbox is not None:
            if not (isinstance(plot_bbox, (tuple, list)) and len(plot_bbox) == 4):
                raise ValueError("plot_bbox must be a tuple (x1, y1, x2, y2)")

            x1, y1, x2, y2 = plot_bbox
            x_min, x_max = min(x1, x2), max(x1, x2)
            y_min, y_max = min(y1, y2), max(y1, y2)

            # boolean masks：兼容 lon/lat 可能为升序或降序
            mask_x = (lon_arr >= x_min) & (lon_arr <= x_max)
            mask_y = (lat_arr >= y_min) & (lat_arr <= y_max)

            idx_x = np.nonzero(mask_x)[0]
            idx_y = np.nonzero(mask_y)[0]

            if idx_x.size == 0 or idx_y.size == 0:
                raise ValueError(f"plot_bbox {plot_bbox} is outside lon/lat range")

            # 使用 np.ix_ 保留原顺序（索引按升序排列），这会截取所需的子网格
            arr_subset = arr[np.ix_(idx_x, idx_y)]
            lon_subset = lon_arr[idx_x]
            lat_subset = lat_arr[idx_y]
            arr = arr_subset  # 替换为裁切后数组
        else:
            lon_subset = lon_arr
            lat_subset = lat_arr

        # 转置以满足 imshow(rows=y, cols=x) 的要求
        data_to_plot = arr.T  # shape -> (len(lat_subset), len(lon_subset))

        # extent 的顺序：(x0, x1, y0, y1) —— 可为升序或降序
        extent = (lon_subset[0], lon_subset[-1], lat_subset[0], lat_subset[-1])

        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(data_to_plot, extent=extent, origin='lower', aspect='auto',
                    vmin=vmin, vmax=vmax, cmap='viridis')
        ax.set_xlabel('lon')
        ax.set_ylabel('lat')
        ax.set_title(f"{time_str} {data_type} {var_name}")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(var_name)

        # overlay points (black 'x', no annotation)
        if points is not None:
            for (lon_p, lat_p) in points:
                # 如果做了裁切，只绘制落在裁切范围内的点
                if plot_bbox is not None:
                    if not (x_min <= lon_p <= x_max and y_min <= lat_p <= y_max):
                        continue
                ax.plot(lon_p, lat_p, marker='x', markersize=4, markeredgewidth=1, color='k')

        plt.tight_layout()
        # 确保保存目录存在
        os.makedirs(folder, exist_ok=True)
        fig.savefig(path, dpi=150)
        plt.close(fig)

        return path


# ------------------ Usage examples ------------------
# Minimal usage (adapt into your transpond_class callback):
#
# saver = NCWriter('out_nc', lon2d, lat2d)
# plotter = FieldPlotter('out_png', lon2d, lat2d)
#
# # for 3D u (nx, nz, ny) save surface level 0:
# saver.save_2d('ignored_time', 'u', u[:,0,:], data_type='original', units='m/s')
# saver.save_2d('ignored_time', 'u', u_da[:,0,:], data_type='da', units='m/s')
# saver.save_2d('ignored_time', 'u', (u_da[:,0,:]-u[:,0,:]), data_type='increment', units='m/s')
#
# # plotting (points are list of (lon, lat))
# plotter.plot_and_save('ignored_time', 'u', u[:,0,:], data_type='original', points=[(105.5, -2.3)])
#
# Files:
# - NetCDF: out_nc/original.nc, out_nc/da.nc, out_nc/increment.nc (variables inside named by var_name)
# - PNG: out_png/<data_type>/<var_name>/<time>_<data_type>_<var_name>.png

if __name__ == "__main__":
    # quick local test
    nx, ny = 30, 40
    lon2d = np.tile(np.linspace(100, 120, nx).reshape(nx, 1), (1, ny))
    lat2d = np.tile(np.linspace(-10, 10, ny).reshape(1, ny), (nx, 1))

    saver = NCWriter('out_nc', lon2d, lat2d)
    plotter = FieldPlotter('out_png', lon2d, lat2d)

    data = np.sin(np.linspace(0, 2 * np.pi, nx)).reshape(nx, 1) * np.cos(np.linspace(0, 2 * np.pi, ny)).reshape(1, ny)
    tstr = '20250608_120000'

    saver.save_2d(tstr, 'u', data, data_type='original', units='m/s')
    saver.save_2d(tstr, 'u', data * 0.9, data_type='da', units='m/s')
    saver.save_2d(tstr, 'u', (data * 0.9 - data), data_type='increment', units='m/s')

    pts = [(105.5, -2.3), (118.2, 5.0)]
    plotter.plot_and_save(tstr, 'u', data, data_type='original', points=pts)

    print('Example outputs written to out_nc/ and out_png/')