from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from utils.paths import DATA_DIR, ensure_project_dirs


TIME_CANDIDATES = ["time", "TIME", "Time", "time_counter", "t", "L", "T", "valid_time"]
LAT_CANDIDATES = ["lat", "latitude", "LAT", "Latitude", "j", "J", "yt_ocean", "y"]
LON_CANDIDATES = ["lon", "longitude", "LON", "Longitude", "i", "I", "xt_ocean", "x"]


def _find_time_name(ds: xr.Dataset) -> str:
    for name in TIME_CANDIDATES:
        if name in ds.coords or name in ds.variables:
            return name
    for name, coord in ds.coords.items():
        if np.issubdtype(np.asarray(coord.values).dtype, np.datetime64):
            return name
    raise ValueError("No time coordinate found.")


def _standardize_time(ds: xr.Dataset) -> xr.Dataset:
    name = _find_time_name(ds)
    if name != "time":
        ds = ds.rename({name: "time"})
    try:
        times = pd.to_datetime(ds["time"].values)
    except Exception:
        times = pd.to_datetime(np.asarray(ds["time"].values, dtype="datetime64[ns]"), errors="coerce")
    if times.hasnans:
        raise ValueError("Observation time coordinate contains invalid values.")
    return ds.assign_coords(time=pd.DatetimeIndex(times))


def _open_many(pattern: str) -> xr.Dataset:
    files = sorted(glob.glob(str(pattern)))
    if not files:
        raise FileNotFoundError(f"No observation files matched: {pattern}")
    datasets = [_standardize_time(xr.open_dataset(path, decode_times=True, mask_and_scale=True)) for path in files]
    return xr.concat(datasets, dim="time", data_vars="minimal", coords="minimal", compat="override")


def _find_var(ds: xr.Dataset, requested: str) -> str:
    if requested in ds.data_vars:
        return requested
    requested_norm = requested.lower().replace("_", "")
    for name in ds.data_vars:
        if name.lower().replace("_", "") == requested_norm:
            return name
    if len(ds.data_vars) == 1:
        return next(iter(ds.data_vars))
    raise KeyError(f"Variable {requested!r} not found. Available: {list(ds.data_vars)}")


def _find_dim_by_name(dims: list[str], candidates: list[str]) -> str | None:
    lower_map = {dim.lower(): dim for dim in dims}
    for item in candidates:
        if item in dims:
            return item
        if item.lower() in lower_map:
            return lower_map[item.lower()]
    return None


def _to_time_yx(ds: xr.Dataset, var_name: str) -> torch.Tensor:
    actual = _find_var(ds, var_name)
    da = ds[actual]
    dims = list(da.dims)
    if "time" not in dims:
        da = da.expand_dims(time=ds["time"])
        dims = list(da.dims)
    non_time = [dim for dim in dims if dim != "time"]
    if len(non_time) < 2:
        raise ValueError(f"Observation variable {actual} does not have two spatial dims: {dims}")
    lat_dim = _find_dim_by_name(non_time, LAT_CANDIDATES)
    lon_dim = _find_dim_by_name(non_time, LON_CANDIDATES)
    if lat_dim is None or lon_dim is None:
        lat_dim, lon_dim = non_time[-2], non_time[-1]
    extra_dims = [dim for dim in non_time if dim not in {lat_dim, lon_dim}]
    for dim in extra_dims:
        da = da.isel({dim: 0})
    da = da.transpose("time", lat_dim, lon_dim)
    return torch.from_numpy(np.asarray(da.values, dtype=np.float32))


def preprocess_observation(config: dict, output_dir: str | Path | None = None) -> Path:
    ensure_project_dirs()
    output_dir = Path(output_dir or config.get("output_dir", DATA_DIR))
    output_dir.mkdir(parents=True, exist_ok=True)
    atm_ds = _open_many(config["atm_glob"])
    ocn_ds = _open_many(config["ocn_glob"])
    common_times = pd.DatetimeIndex(atm_ds["time"].values).intersection(pd.DatetimeIndex(ocn_ds["time"].values))
    common_times = common_times.sort_values()
    if len(common_times) == 0:
        raise ValueError("Atmosphere and ocean observations have no overlapping times.")
    atm_ds = atm_ds.sel(time=common_times)
    ocn_ds = ocn_ds.sel(time=common_times)
    payload = {
        "data": {
            config.get("atm_output_key", "ps1"): _to_time_yx(atm_ds, config.get("atm_var", "ps1")),
            config.get("ocn_output_key", "sst"): _to_time_yx(ocn_ds, config.get("ocn_var", "sst")),
        },
        "coords": {
            "time": np.asarray([pd.Timestamp(t).isoformat() for t in common_times], dtype=object),
        },
        "shapes": {},
    }
    payload["shapes"] = {key: tuple(value.shape) for key, value in payload["data"].items()}
    out_path = output_dir / config.get("output_name", "merged_observation.pt")
    torch.save(payload, out_path)
    return out_path


def _masked_observation(data: torch.Tensor) -> np.ma.MaskedArray:
    arr = np.asarray(data.detach().cpu().numpy(), dtype=np.float64)
    return np.ma.masked_invalid(arr)


def plot_observation_file(
    observation_path: str | Path,
    output_dir: str | Path,
    max_times: int | None = None,
) -> None:
    payload = torch.load(observation_path, map_location="cpu", weights_only=False)
    data = payload["data"]
    times = payload.get("coords", {}).get("time", [])
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, tensor in data.items():
        tensor = tensor.squeeze()
        if tensor.dim() != 3:
            continue
        count = tensor.shape[0] if max_times is None else min(int(max_times), tensor.shape[0])
        var_dir = output_dir / key
        var_dir.mkdir(parents=True, exist_ok=True)
        for idx in range(count):
            field = _masked_observation(tensor[idx])
            values = field.compressed()
            if values.size == 0:
                continue
            fig, ax = plt.subplots(figsize=(8.0, 4.8))
            vmin = float(np.nanpercentile(values, 0.5))
            vmax = float(np.nanpercentile(values, 99.5))
            cmap = "turbo" if key.lower() == "sst" else "bwr"
            image = ax.imshow(field, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            label = str(times[idx]) if idx < len(times) else f"index={idx}"
            ax.set_title(f"{key} observation {idx} {label}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            fig.colorbar(image, ax=ax, label=key)
            fig.tight_layout()
            fig.savefig(var_dir / f"{idx:04d}_{key}_observation.png", dpi=150)
            plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "preprocess_observation.json"))
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--plot-only", default=None)
    parser.add_argument("--plot-output-dir", default=str(PROJECT_DIR / "plots" / "observation_check"))
    parser.add_argument("--plot-max-times", type=int, default=24)
    args = parser.parse_args()
    if args.plot_only:
        plot_observation_file(args.plot_only, args.plot_output_dir, args.plot_max_times)
        print(f"saved observation plots to {args.plot_output_dir}")
        return
    with open(args.config, "r", encoding="utf-8") as f:
        config = json.load(f)
    print(f"saved {preprocess_observation(config, args.output_dir)}")


if __name__ == "__main__":
    main()
