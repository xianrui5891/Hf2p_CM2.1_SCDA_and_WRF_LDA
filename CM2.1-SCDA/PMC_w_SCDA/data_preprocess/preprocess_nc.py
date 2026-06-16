from __future__ import annotations

import argparse
import gc
import glob
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

import numpy as np
import torch
import xarray as xr

try:
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover - tqdm is optional but recommended.
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else range(kwargs.get("total", 0))


PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from utils.data import RUNTIME_ATM_VARS, RUNTIME_OCN_VARS
from utils.paths import DATA_DIR, ensure_project_dirs


TIME_NAMES = ("time", "TIME", "Time", "time_counter", "t", "L", "T", "valid_time")
VERT_NAMES = ("K", "k", "lev", "level", "z", "Z", "depth", "pfull", "phalf", "st_ocean", "st_ocean_levels", "zt")
DEFAULT_VAR_ALIASES = {
    "atm": {
        "T_SURF": ["T_SURF", "TEMP", "temp", "t_surf", "tsurf"],
        "U_REF": ["U_REF", "UCOMP", "u_ref", "ucomp"],
        "V_REF": ["V_REF", "VCOMP", "v_ref", "vcomp"],
        "PS": ["PS", "ps", "PRES", "pressure"],
        "TEMP": ["TEMP", "T_SURF", "temp", "t_surf"],
        "UCOMP": ["UCOMP", "U_REF", "ucomp", "u_ref"],
        "VCOMP": ["VCOMP", "V_REF", "vcomp", "v_ref"],
    },
    "ocn": {
        "SST": ["SST", "TEMP", "sst", "temp"],
        "U_SURF": ["U_SURF", "U", "u_surf", "u"],
        "V_SURF": ["V_SURF", "V", "v_surf", "v"],
        "ETA_T": ["ETA_T", "eta_t", "SSH", "ssh", "SSZ", "ssz"],
    },
}
_PROGRESS_LOG_DIR: Path | None = None


def _set_progress_log_dir(path: str | Path | None) -> None:
    global _PROGRESS_LOG_DIR
    if path is None:
        _PROGRESS_LOG_DIR = None
        return
    _PROGRESS_LOG_DIR = Path(path)
    _PROGRESS_LOG_DIR.mkdir(parents=True, exist_ok=True)


def _progress_log(name: str, message: str) -> None:
    text = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
    print(text, flush=True)
    if _PROGRESS_LOG_DIR is not None:
        with open(_PROGRESS_LOG_DIR / name, "a", encoding="utf-8") as f:
            f.write(text + "\n")


@dataclass
class FilePair:
    atm_path: Path
    ocn_path: Path
    length: int
    global_start: int = 0


class RunningStats:
    def __init__(self) -> None:
        self.total = 0.0
        self.total_sq = 0.0
        self.count = 0

    def update(self, values: np.ndarray) -> None:
        finite = np.asarray(values, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        self.total += float(finite.sum())
        self.total_sq += float(np.square(finite).sum())
        self.count += int(finite.size)

    def merge(self, other: "RunningStats") -> None:
        self.total += other.total
        self.total_sq += other.total_sq
        self.count += other.count

    def to_raw(self) -> tuple[float, float, int]:
        return self.total, self.total_sq, self.count

    @classmethod
    def from_raw(cls, raw: tuple[float, float, int]) -> "RunningStats":
        item = cls()
        item.total, item.total_sq, item.count = float(raw[0]), float(raw[1]), int(raw[2])
        return item

    def to_dict(self) -> dict[str, float]:
        if self.count == 0:
            return {"mean": 0.0, "std": 1.0}
        mean = self.total / self.count
        variance = max(self.total_sq / self.count - mean * mean, 1e-16)
        return {"mean": float(mean), "std": float(math.sqrt(variance))}


def _list_files(pattern: str) -> list[Path]:
    files = [Path(path) for path in sorted(glob.glob(str(pattern)))]
    if not files:
        raise FileNotFoundError(f"No NetCDF files matched: {pattern}")
    return files


def _get_mpi_comm():
    try:
        from mpi4py import MPI
    except Exception:
        return None
    return MPI.COMM_WORLD


def _open_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(
        path,
        decode_times=True,
        decode_timedelta=True,
        mask_and_scale=True,
        cache=False,
    )


def _find_time_dim(da: xr.DataArray) -> str | None:
    for name in TIME_NAMES:
        if name in da.dims:
            return name
    return None


def _find_spatial_dims(da: xr.DataArray, time_dim: str | None, vert_dim: str | None) -> tuple[str, str]:
    dims = [dim for dim in da.dims if dim not in {time_dim, vert_dim}]
    if len(dims) < 2:
        raise ValueError(f"Cannot find two spatial dims for {da.name}: {da.dims}")
    return dims[-2], dims[-1]


def _find_vert_dim(da: xr.DataArray, time_dim: str | None) -> str | None:
    candidates = [dim for dim in da.dims if dim != time_dim]
    if len(candidates) <= 2:
        return None
    for name in VERT_NAMES:
        if name in da.dims:
            return name
    spatial_names = {"J", "I", "j", "i", "lat", "lon", "latitude", "longitude", "y", "x", "yt_ocean", "xt_ocean"}
    for dim in candidates:
        if dim not in spatial_names:
            return dim
    return candidates[0]


def _resolve_level(selection: int, size: int) -> int:
    if selection == -1:
        return 0
    if selection == -2:
        return size - 1
    if selection < 0:
        return max(0, size + selection)
    return min(selection, size - 1)


def _normalize_name(name: str) -> str:
    return name.lower().replace("_", "").replace("-", "")


def _resolve_var_name(ds: xr.Dataset, var_name: str, domain: str, var_aliases: dict) -> str:
    if var_name in ds.data_vars:
        return var_name
    var_lower = var_name.lower()
    for name in ds.data_vars:
        if name.lower() == var_lower:
            return name
    var_norm = _normalize_name(var_name)
    for name in ds.data_vars:
        if _normalize_name(name) == var_norm:
            return name

    aliases = []
    aliases.extend(DEFAULT_VAR_ALIASES.get(domain, {}).get(var_name, []))
    aliases.extend(DEFAULT_VAR_ALIASES.get(domain, {}).get(var_name.upper(), []))
    aliases.extend(var_aliases.get(domain, {}).get(var_name, []))
    aliases.extend(var_aliases.get(domain, {}).get(var_name.upper(), []))
    for alias in aliases:
        if alias in ds.data_vars:
            return alias
        alias_norm = _normalize_name(str(alias))
        for name in ds.data_vars:
            if _normalize_name(name) == alias_norm:
                return name

    available = ", ".join(list(ds.data_vars)[:40])
    raise KeyError(f"{var_name!r} not found in {domain} dataset. Available variables: {available}")


def _prepare_var(ds: xr.Dataset, var_name: str, domain: str, vert_select: dict, var_aliases: dict) -> xr.DataArray:
    actual_name = _resolve_var_name(ds, var_name, domain, var_aliases)
    da = ds[actual_name]
    time_dim = _find_time_dim(da)
    vert_dim = _find_vert_dim(da, time_dim)
    if vert_dim is not None:
        key = f"{domain}::{var_name}"
        actual_key = f"{domain}::{actual_name}"
        level = _resolve_level(
            int(vert_select.get(key, vert_select.get(actual_key, vert_select.get(var_name, vert_select.get(actual_name, -1))))),
            da.sizes[vert_dim],
        )
        da = da.isel({vert_dim: level})
    y_dim, x_dim = _find_spatial_dims(da, time_dim, None)
    if time_dim is None:
        da = da.expand_dims(time=[0])
        time_dim = "time"
    return da.transpose(time_dim, y_dim, x_dim)


def _domain_length(ds: xr.Dataset, vars_: Sequence[str], domain: str, vert_select: dict, var_aliases: dict) -> int:
    lengths = []
    for name in vars_:
        da = _prepare_var(ds, name, domain, vert_select, var_aliases)
        lengths.append(int(da.sizes[_find_time_dim(da) or "time"]))
    return min(lengths)


def _dataset_time_values(ds: xr.Dataset, length: int) -> list[str | int | float]:
    time_name = next((name for name in TIME_NAMES if name in ds.coords or name in ds.variables), None)
    if time_name is None:
        return list(range(length))
    values = np.asarray(ds[time_name].values)[:length]
    out: list[str | int | float] = []
    for value in values:
        if np.issubdtype(np.asarray(value).dtype, np.datetime64):
            out.append(str(np.datetime_as_string(value, unit="s")))
        else:
            try:
                out.append(float(value))
            except Exception:
                out.append(str(value))
    return out


def _read_var_chunk(
    ds: xr.Dataset,
    var_name: str,
    domain: str,
    vert_select: dict,
    var_aliases: dict,
    start: int,
    stop: int,
) -> np.ndarray:
    da = _prepare_var(ds, var_name, domain, vert_select, var_aliases)
    time_dim = _find_time_dim(da) or "time"
    chunk = da.isel({time_dim: slice(start, stop)})
    return np.asarray(chunk.values, dtype=np.float32)


def _empty_stats(vars_: Sequence[str], domain: str) -> dict[str, RunningStats]:
    return {f"{domain}::{name}": RunningStats() for name in vars_}


def _scan_and_accumulate_stats(
    atm_files: Sequence[Path],
    ocn_files: Sequence[Path],
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    time_chunk: int,
) -> tuple[list[FilePair], dict[str, dict[str, float]]]:
    """Build the ordered file plan and z-score stats in one streaming pass."""
    stats = {**_empty_stats(atm_output_vars, "atm"), **_empty_stats(ocn_output_vars, "ocn")}
    pairs: list[FilePair] = []
    frame_progress = tqdm(desc="pass 1/2 scan+z-score", unit="frame")
    file_progress = tqdm(zip(atm_files, ocn_files), total=min(len(atm_files), len(ocn_files)), desc="source files", unit="file")
    try:
        for atm_path, ocn_path in file_progress:
            with _open_dataset(atm_path) as atm_ds, _open_dataset(ocn_path) as ocn_ds:
                length = min(
                    _domain_length(atm_ds, atm_vars, "atm", vert_select, var_aliases),
                    _domain_length(ocn_ds, ocn_vars, "ocn", vert_select, var_aliases),
                )
                if length <= 0:
                    continue
                pairs.append(FilePair(atm_path=atm_path, ocn_path=ocn_path, length=length))

                for start in range(0, length, time_chunk):
                    stop = min(start + time_chunk, length)
                    for raw_name, out_name in zip(atm_vars, atm_output_vars):
                        chunk = _read_var_chunk(atm_ds, raw_name, "atm", vert_select, var_aliases, start, stop)
                        stats[f"atm::{out_name}"].update(chunk)
                        del chunk
                    for raw_name, out_name in zip(ocn_vars, ocn_output_vars):
                        chunk = _read_var_chunk(ocn_ds, raw_name, "ocn", vert_select, var_aliases, start, stop)
                        stats[f"ocn::{out_name}"].update(chunk)
                        del chunk
                    frame_progress.update(stop - start)
            gc.collect()
    finally:
        frame_progress.close()
        file_progress.close()
    return pairs, {key: value.to_dict() for key, value in stats.items()}


def _chunk_sequence(items: Sequence[FilePair], chunks: int) -> list[list[FilePair]]:
    chunks = max(1, int(chunks))
    out: list[list[FilePair]] = [[] for _ in range(chunks)]
    loads = [0 for _ in range(chunks)]
    for item in items:
        idx = int(np.argmin(loads))
        out[idx].append(item)
        loads[idx] += int(item.length)
    return [chunk for chunk in out if chunk]


def _chunk_indexed_pairs(
    atm_files: Sequence[Path],
    ocn_files: Sequence[Path],
    chunks: int,
) -> list[list[tuple[int, Path, Path]]]:
    chunks = max(1, int(chunks))
    out: list[list[tuple[int, Path, Path]]] = [[] for _ in range(chunks)]
    for idx, (atm_path, ocn_path) in enumerate(zip(atm_files, ocn_files)):
        out[idx % chunks].append((idx, Path(atm_path), Path(ocn_path)))
    return [chunk for chunk in out if chunk]


def _finalize_plan_entries(
    entries: Sequence[tuple[int, str, str, int, list[str | int | float]]],
) -> tuple[list[FilePair], list[str | int | float]]:
    pairs: list[FilePair] = []
    time_values: list[str | int | float] = []
    global_start = 0
    for _index, atm_path, ocn_path, length, item_times in sorted(entries, key=lambda item: item[0]):
        if int(length) <= 0:
            continue
        pairs.append(
            FilePair(
                atm_path=Path(atm_path),
                ocn_path=Path(ocn_path),
                length=int(length),
                global_start=global_start,
            )
        )
        time_values.extend(item_times)
        global_start += int(length)
    return pairs, time_values


def _stats_and_plan_worker(
    worker_id: int,
    indexed_pairs: Sequence[tuple[int, Path, Path]],
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    time_chunk: int,
) -> dict[str, object]:
    stats = {**_empty_stats(atm_output_vars, "atm"), **_empty_stats(ocn_output_vars, "ocn")}
    plan_entries: list[tuple[int, str, str, int, list[str | int | float]]] = []
    done_files = 0
    done_frames = 0
    total_files = len(indexed_pairs)
    log_name = f"preprocess_nc_rank{worker_id:02d}.log"
    _progress_log(log_name, f"[stats+plan rank={worker_id}] start files={total_files}")
    for file_index, atm_path, ocn_path in indexed_pairs:
        with _open_dataset(atm_path) as atm_ds, _open_dataset(ocn_path) as ocn_ds:
            length = min(
                _domain_length(atm_ds, atm_vars, "atm", vert_select, var_aliases),
                _domain_length(ocn_ds, ocn_vars, "ocn", vert_select, var_aliases),
            )
            if length > 0:
                item_times = _dataset_time_values(atm_ds, length)
                plan_entries.append((int(file_index), str(atm_path), str(ocn_path), int(length), item_times))
                for start in range(0, length, time_chunk):
                    stop = min(start + time_chunk, length)
                    for raw_name, out_name in zip(atm_vars, atm_output_vars):
                        chunk = _read_var_chunk(atm_ds, raw_name, "atm", vert_select, var_aliases, start, stop)
                        stats[f"atm::{out_name}"].update(chunk)
                        del chunk
                    for raw_name, out_name in zip(ocn_vars, ocn_output_vars):
                        chunk = _read_var_chunk(ocn_ds, raw_name, "ocn", vert_select, var_aliases, start, stop)
                        stats[f"ocn::{out_name}"].update(chunk)
                        del chunk
                    done_frames += stop - start
                    if done_frames % 500 == 0:
                        _progress_log(log_name, f"[stats+plan rank={worker_id}] frames={done_frames}")
        done_files += 1
        if done_files == total_files or done_files % 10 == 0:
            _progress_log(log_name, f"[stats+plan rank={worker_id}] files={done_files}/{total_files} frames={done_frames}")
        gc.collect()
    _progress_log(log_name, f"[stats+plan rank={worker_id}] done files={done_files} frames={done_frames}")
    return {
        "stats": {key: value.to_raw() for key, value in stats.items()},
        "plan_entries": plan_entries,
    }


def _merge_raw_stats(raw_items: Sequence[dict[str, tuple[float, float, int]]]) -> dict[str, dict[str, float]]:
    merged: dict[str, RunningStats] = {}
    for raw in raw_items:
        for key, value in raw.items():
            merged.setdefault(key, RunningStats()).merge(RunningStats.from_raw(value))
    return {key: value.to_dict() for key, value in merged.items()}


def _normalize_chunk(
    data: np.ndarray,
    vars_: Sequence[str],
    domain: str,
    normalization_stats: dict[str, dict[str, float]],
) -> np.ndarray:
    out = np.asarray(data, dtype=np.float32).copy()
    for idx, name in enumerate(vars_):
        stat = normalization_stats[f"{domain}::{name}"]
        out[:, idx] = (out[:, idx] - float(stat["mean"])) / max(float(stat["std"]), 1e-8)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def _read_domain_chunk(
    ds: xr.Dataset,
    raw_vars: Sequence[str],
    domain: str,
    vert_select: dict,
    var_aliases: dict,
    start: int,
    stop: int,
) -> np.ndarray:
    arrays = [_read_var_chunk(ds, name, domain, vert_select, var_aliases, start, stop) for name in raw_vars]
    return np.stack(arrays, axis=1)


def _cleanup_old_frames(output_dir: Path) -> None:
    for path in output_dir.glob("frame_*.pt"):
        path.unlink()


def _run_parallel_stats_and_plan(
    atm_files: Sequence[Path],
    ocn_files: Sequence[Path],
    workers: int,
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    time_chunk: int,
) -> tuple[list[FilePair], list[str | int | float], dict[str, dict[str, float]]]:
    chunks = _chunk_indexed_pairs(atm_files, ocn_files, workers)
    results: list[dict[str, object]] = []
    progress = tqdm(total=len(chunks), desc="pass 1/2 stats+plan workers", unit="worker")
    try:
        with ProcessPoolExecutor(max_workers=len(chunks)) as executor:
            futures = [
                executor.submit(
                    _stats_and_plan_worker,
                    idx,
                    chunk,
                    atm_vars,
                    atm_output_vars,
                    ocn_vars,
                    ocn_output_vars,
                    vert_select,
                    var_aliases,
                    time_chunk,
                )
                for idx, chunk in enumerate(chunks)
            ]
            for future in as_completed(futures):
                results.append(future.result())
                progress.update(1)
    finally:
        progress.close()
    raw_stats = [item["stats"] for item in results]
    plan_entries = [
        entry
        for item in results
        for entry in item["plan_entries"]
    ]
    plan, time_values = _finalize_plan_entries(plan_entries)
    return plan, time_values, _merge_raw_stats(raw_stats)


def _run_parallel_write_frames(
    plan: Sequence[FilePair],
    output_dir: Path,
    workers: int,
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    normalization_stats: dict[str, dict[str, float]],
    time_chunk: int,
) -> tuple[list[str], list[int], list[int]]:
    chunks = _chunk_sequence(plan, workers)
    total_frames = sum(item.length for item in plan)
    progress = tqdm(total=len(chunks), desc="pass 3/3 write workers", unit="worker")
    results: list[dict[str, object]] = []
    try:
        with ProcessPoolExecutor(max_workers=len(chunks)) as executor:
            futures = [
                executor.submit(
                    _write_frames_worker,
                    idx,
                    chunk,
                    output_dir,
                    atm_vars,
                    atm_output_vars,
                    ocn_vars,
                    ocn_output_vars,
                    vert_select,
                    var_aliases,
                    normalization_stats,
                    time_chunk,
                )
                for idx, chunk in enumerate(chunks)
            ]
            for future in as_completed(futures):
                results.append(future.result())
                progress.update(1)
    finally:
        progress.close()
    written = sum(int(item.get("written", 0)) for item in results)
    if written != total_frames:
        raise RuntimeError(f"Parallel preprocessing wrote {written} frames, expected {total_frames}.")
    first_atm = next((item.get("atm_shape") for item in results if item.get("atm_shape") is not None), None)
    first_ocn = next((item.get("ocn_shape") for item in results if item.get("ocn_shape") is not None), None)
    atm_shape = list(first_atm) if first_atm is not None else [0, len(atm_output_vars), 0, 0]
    ocn_shape = list(first_ocn) if first_ocn is not None else [0, len(ocn_output_vars), 0, 0]
    atm_shape[0] = total_frames
    ocn_shape[0] = total_frames
    frame_files = [f"frame_{idx:08d}.pt" for idx in range(total_frames)]
    return frame_files, atm_shape, ocn_shape


def _preprocess_mpi(
    config: dict,
    output_dir: Path,
    comm,
) -> tuple[Path | None, Path]:
    rank = int(comm.Get_rank())
    world_size = int(comm.Get_size())

    if rank == 0:
        ensure_project_dirs()
        output_dir.mkdir(parents=True, exist_ok=True)
        _set_progress_log_dir(config.get("progress_log_dir", output_dir.parent / "logs"))
        if bool(config.get("cleanup_existing_frames", True)):
            _cleanup_old_frames(output_dir)
        atm_files = _list_files(config["atm_glob"])
        ocn_files = _list_files(config["ocn_glob"])
        if len(atm_files) != len(ocn_files):
            print(f"[WARN] ATM file count {len(atm_files)} != OCN file count {len(ocn_files)}. Pairing by sorted order.", flush=True)
    else:
        atm_files = []
        ocn_files = []

    atm_files = comm.bcast(atm_files, root=0)
    ocn_files = comm.bcast(ocn_files, root=0)
    output_dir = Path(comm.bcast(str(output_dir), root=0))
    progress_log_dir = comm.bcast(str(config.get("progress_log_dir", output_dir.parent / "logs")), root=0)
    _set_progress_log_dir(progress_log_dir)

    atm_vars = list(config.get("atm_vars", RUNTIME_ATM_VARS))
    atm_output_vars = list(config.get("atm_output_vars", atm_vars))
    ocn_vars = list(config.get("ocn_vars", RUNTIME_OCN_VARS))
    ocn_output_vars = list(config.get("ocn_output_vars", RUNTIME_OCN_VARS))
    _validate_config_vars(atm_vars, atm_output_vars, "atm")
    _validate_config_vars(ocn_vars, ocn_output_vars, "ocn")
    vert_select = dict(config.get("vert_select", {}))
    var_aliases = dict(config.get("var_aliases", {}))
    time_chunk = max(1, int(config.get("time_chunk", 1)))

    indexed_chunks = _chunk_indexed_pairs(atm_files, ocn_files, world_size)
    local_indexed_pairs = indexed_chunks[rank] if rank < len(indexed_chunks) else []
    local_stats_plan = _stats_and_plan_worker(
        rank,
        local_indexed_pairs,
        atm_vars,
        atm_output_vars,
        ocn_vars,
        ocn_output_vars,
        vert_select,
        var_aliases,
        time_chunk,
    )
    gathered_stats_plan = comm.gather(local_stats_plan, root=0)
    if rank == 0:
        normalization_stats = _merge_raw_stats([item["stats"] for item in gathered_stats_plan])
        plan_entries = [
            entry
            for item in gathered_stats_plan
            for entry in item["plan_entries"]
        ]
        plan, time_values = _finalize_plan_entries(plan_entries)
        if not plan:
            raise RuntimeError("No aligned frames found in the configured NetCDF files.")
        _progress_log(
            "preprocess_nc_rank00.log",
            f"[preprocess mpi] stats+plan files={len(plan)} frames={sum(item.length for item in plan)} ranks={world_size}",
        )
        out_stats = output_dir / "normalization_stats.json"
        with open(out_stats, "w", encoding="utf-8") as f:
            json.dump(normalization_stats, f, indent=2)
    else:
        normalization_stats = None
        plan = None
        time_values = None
    plan = comm.bcast(plan, root=0)
    time_values = comm.bcast(time_values, root=0)
    normalization_stats = comm.bcast(normalization_stats, root=0)
    comm.Barrier()

    chunks = _chunk_sequence(plan, world_size)
    local_plan = chunks[rank] if rank < len(chunks) else []

    local_result = _write_frames_worker(
        rank,
        local_plan,
        output_dir,
        atm_vars,
        atm_output_vars,
        ocn_vars,
        ocn_output_vars,
        vert_select,
        var_aliases,
        normalization_stats,
        time_chunk,
    )
    gathered_results = comm.gather(local_result, root=0)
    if rank != 0:
        return None, output_dir / "normalization_metadata.json"

    total_frames = sum(item.length for item in plan)
    written = sum(int(item.get("written", 0)) for item in gathered_results)
    if written != total_frames:
        raise RuntimeError(f"MPI preprocessing wrote {written} frames, expected {total_frames}.")
    first_atm = next((item.get("atm_shape") for item in gathered_results if item.get("atm_shape") is not None), None)
    first_ocn = next((item.get("ocn_shape") for item in gathered_results if item.get("ocn_shape") is not None), None)
    atm_shape = list(first_atm) if first_atm is not None else [0, len(atm_output_vars), 0, 0]
    ocn_shape = list(first_ocn) if first_ocn is not None else [0, len(ocn_output_vars), 0, 0]
    atm_shape[0] = total_frames
    ocn_shape[0] = total_frames
    frame_files = [f"frame_{idx:08d}.pt" for idx in range(total_frames)]
    out_meta = output_dir / "normalization_metadata.json"
    out_stats = output_dir / "normalization_stats.json"
    metadata = {
        "atm_vars": atm_output_vars,
        "raw_atm_vars": atm_vars,
        "ocn_vars": ocn_output_vars,
        "raw_ocn_vars": ocn_vars,
        "var_aliases": var_aliases,
        "normalization": "palm_gda_zscore",
        "normalization_stats_file": out_stats.name,
        "normalization_stats_path": str(out_stats),
        "handle_nan": True,
        "nan_fill_value": 0.0,
        "atm_shape": atm_shape,
        "ocn_shape": ocn_shape,
        "vert_mode": "fold_z",
        "streaming": True,
        "parallel_backend": "mpi4py",
        "preprocess_ranks": world_size,
        "time_chunk": time_chunk,
        "num_source_pairs": len(plan),
        "num_frames": total_frames,
        "frame_files": frame_files,
        "time_values": list(time_values),
        "aggregate_file": None,
    }
    with open(out_meta, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    with open(output_dir / "_frame_manifest.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    return None, out_meta


def _write_frames(
    plan: Sequence[FilePair],
    output_dir: Path,
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    normalization_stats: dict[str, dict[str, float]],
    time_chunk: int,
) -> tuple[list[str], list[str | int | float], list[int], list[int]]:
    frame_files: list[str] = []
    time_values: list[str | int | float] = []
    atm_shape: list[int] | None = None
    ocn_shape: list[int] | None = None
    global_idx = 0
    total_frames = sum(item.length for item in plan)
    progress = tqdm(total=total_frames, desc="pass 2/2 write frames", unit="frame")
    try:
        for item in plan:
            with _open_dataset(item.atm_path) as atm_ds, _open_dataset(item.ocn_path) as ocn_ds:
                file_time_values = _dataset_time_values(atm_ds, item.length)
                for start in range(0, item.length, time_chunk):
                    stop = min(start + time_chunk, item.length)
                    atm = _read_domain_chunk(atm_ds, atm_vars, "atm", vert_select, var_aliases, start, stop)
                    ocn = _read_domain_chunk(ocn_ds, ocn_vars, "ocn", vert_select, var_aliases, start, stop)
                    atm_mask = np.isfinite(atm)
                    ocn_mask = np.isfinite(ocn)
                    atm = _normalize_chunk(atm, atm_output_vars, "atm", normalization_stats)
                    ocn = _normalize_chunk(ocn, ocn_output_vars, "ocn", normalization_stats)
                    if atm_shape is None:
                        atm_shape = [0, *list(atm.shape[1:])]
                    if ocn_shape is None:
                        ocn_shape = [0, *list(ocn.shape[1:])]
                    for local_idx in range(stop - start):
                        filename = f"frame_{global_idx:08d}.pt"
                        torch.save(
                            {
                                "atm": torch.from_numpy(atm[local_idx]).contiguous(),
                                "ocn": torch.from_numpy(ocn[local_idx]).contiguous(),
                                "atm_mask": torch.from_numpy(atm_mask[local_idx]).contiguous(),
                                "ocn_mask": torch.from_numpy(ocn_mask[local_idx]).contiguous(),
                                "index": global_idx,
                                "time_index": global_idx,
                                "time": file_time_values[start + local_idx] if start + local_idx < len(file_time_values) else global_idx,
                                "source_atm_file": item.atm_path.name,
                                "source_ocn_file": item.ocn_path.name,
                            },
                            output_dir / filename,
                        )
                        frame_files.append(filename)
                        time_values.append(file_time_values[start + local_idx] if start + local_idx < len(file_time_values) else global_idx)
                        global_idx += 1
                        progress.update(1)
                    del atm, ocn, atm_mask, ocn_mask
                del file_time_values
            gc.collect()
    finally:
        progress.close()
    if atm_shape is None:
        atm_shape = [0, len(atm_output_vars), 0, 0]
    if ocn_shape is None:
        ocn_shape = [0, len(ocn_output_vars), 0, 0]
    atm_shape[0] = len(frame_files)
    ocn_shape[0] = len(frame_files)
    return frame_files, time_values, atm_shape, ocn_shape


def _write_frames_worker(
    worker_id: int,
    plan: Sequence[FilePair],
    output_dir: str | Path,
    atm_vars: Sequence[str],
    atm_output_vars: Sequence[str],
    ocn_vars: Sequence[str],
    ocn_output_vars: Sequence[str],
    vert_select: dict,
    var_aliases: dict,
    normalization_stats: dict[str, dict[str, float]],
    time_chunk: int,
) -> dict[str, object]:
    output_dir = Path(output_dir)
    atm_shape: list[int] | None = None
    ocn_shape: list[int] | None = None
    written = 0
    total_frames = sum(item.length for item in plan)
    log_name = f"preprocess_nc_rank{worker_id:02d}.log"
    _progress_log(log_name, f"[write rank={worker_id}] start files={len(plan)} frames={total_frames}")
    for item in plan:
        with _open_dataset(item.atm_path) as atm_ds, _open_dataset(item.ocn_path) as ocn_ds:
            file_time_values = _dataset_time_values(atm_ds, item.length)
            for start in range(0, item.length, time_chunk):
                stop = min(start + time_chunk, item.length)
                atm = _read_domain_chunk(atm_ds, atm_vars, "atm", vert_select, var_aliases, start, stop)
                ocn = _read_domain_chunk(ocn_ds, ocn_vars, "ocn", vert_select, var_aliases, start, stop)
                atm_mask = np.isfinite(atm)
                ocn_mask = np.isfinite(ocn)
                atm = _normalize_chunk(atm, atm_output_vars, "atm", normalization_stats)
                ocn = _normalize_chunk(ocn, ocn_output_vars, "ocn", normalization_stats)
                if atm_shape is None:
                    atm_shape = [0, *list(atm.shape[1:])]
                if ocn_shape is None:
                    ocn_shape = [0, *list(ocn.shape[1:])]
                for local_idx in range(stop - start):
                    global_idx = item.global_start + start + local_idx
                    filename = f"frame_{global_idx:08d}.pt"
                    torch.save(
                        {
                            "atm": torch.from_numpy(atm[local_idx]).contiguous(),
                            "ocn": torch.from_numpy(ocn[local_idx]).contiguous(),
                            "atm_mask": torch.from_numpy(atm_mask[local_idx]).contiguous(),
                            "ocn_mask": torch.from_numpy(ocn_mask[local_idx]).contiguous(),
                            "index": global_idx,
                            "time_index": global_idx,
                            "time": file_time_values[start + local_idx] if start + local_idx < len(file_time_values) else global_idx,
                            "source_atm_file": item.atm_path.name,
                            "source_ocn_file": item.ocn_path.name,
                        },
                        output_dir / filename,
                    )
                    written += 1
                    if written == total_frames or written % 500 == 0:
                        _progress_log(log_name, f"[write rank={worker_id}] frames={written}/{total_frames}")
                del atm, ocn, atm_mask, ocn_mask
            del file_time_values
        gc.collect()
    _progress_log(log_name, f"[write rank={worker_id}] done files={len(plan)} frames={written}")
    return {"written": written, "atm_shape": atm_shape, "ocn_shape": ocn_shape}


def _validate_config_vars(raw_vars: Sequence[str], output_vars: Sequence[str], domain: str) -> None:
    if len(raw_vars) != len(output_vars):
        raise ValueError(f"{domain} raw vars and output vars must have the same length.")


def preprocess(config: dict, output_dir: str | Path = DATA_DIR) -> tuple[Path | None, Path]:
    _set_progress_log_dir(config.get("progress_log_dir"))
    mpi_comm = _get_mpi_comm()
    parallel_backend = str(config.get("parallel_backend", "auto")).lower()
    if mpi_comm is not None and int(mpi_comm.Get_size()) > 1 and parallel_backend in {"auto", "mpi", "mpi4py"}:
        return _preprocess_mpi(config, Path(output_dir), mpi_comm)

    ensure_project_dirs()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if bool(config.get("cleanup_existing_frames", True)):
        _cleanup_old_frames(output_dir)

    atm_files = _list_files(config["atm_glob"])
    ocn_files = _list_files(config["ocn_glob"])
    if len(atm_files) != len(ocn_files):
        print(f"[WARN] ATM file count {len(atm_files)} != OCN file count {len(ocn_files)}. Pairing by sorted order.")

    atm_vars = list(config.get("atm_vars", RUNTIME_ATM_VARS))
    atm_output_vars = list(config.get("atm_output_vars", atm_vars))
    ocn_vars = list(config.get("ocn_vars", RUNTIME_OCN_VARS))
    ocn_output_vars = list(config.get("ocn_output_vars", RUNTIME_OCN_VARS))
    _validate_config_vars(atm_vars, atm_output_vars, "atm")
    _validate_config_vars(ocn_vars, ocn_output_vars, "ocn")
    vert_select = dict(config.get("vert_select", {}))
    var_aliases = dict(config.get("var_aliases", {}))
    time_chunk = max(1, int(config.get("time_chunk", 1)))
    preprocess_workers = max(1, int(config.get("preprocess_workers", config.get("num_workers", 1))))
    parallel_backend = str(config.get("parallel_backend", "auto")).lower()

    if preprocess_workers > 1 and parallel_backend in {"auto", "process", "multiprocessing", "local"}:
        plan, time_values, normalization_stats = _run_parallel_stats_and_plan(
            atm_files,
            ocn_files,
            preprocess_workers,
            atm_vars,
            atm_output_vars,
            ocn_vars,
            ocn_output_vars,
            vert_select,
            var_aliases,
            time_chunk=time_chunk,
        )
    else:
        plan, normalization_stats = _scan_and_accumulate_stats(
            atm_files,
            ocn_files,
            atm_vars,
            atm_output_vars,
            ocn_vars,
            ocn_output_vars,
            vert_select,
            var_aliases,
            time_chunk=time_chunk,
        )
        time_values = None
    if not plan:
        raise RuntimeError("No aligned frames found in the configured NetCDF files.")
    out_meta = output_dir / "normalization_metadata.json"
    out_stats = output_dir / "normalization_stats.json"
    with open(out_stats, "w", encoding="utf-8") as f:
        json.dump(normalization_stats, f, indent=2)

    if preprocess_workers > 1 and parallel_backend in {"auto", "process", "multiprocessing", "local"}:
        frame_files, atm_shape, ocn_shape = _run_parallel_write_frames(
            plan,
            output_dir,
            preprocess_workers,
            atm_vars,
            atm_output_vars,
            ocn_vars,
            ocn_output_vars,
            vert_select,
            var_aliases,
            normalization_stats,
            time_chunk=time_chunk,
        )
        if time_values is None:
            raise RuntimeError("Parallel preprocessing did not produce time_values.")
    else:
        frame_files, time_values, atm_shape, ocn_shape = _write_frames(
            plan,
            output_dir,
            atm_vars,
            atm_output_vars,
            ocn_vars,
            ocn_output_vars,
            vert_select,
            var_aliases,
            normalization_stats,
            time_chunk=time_chunk,
        )

    metadata = {
        "atm_vars": atm_output_vars,
        "raw_atm_vars": atm_vars,
        "ocn_vars": ocn_output_vars,
        "raw_ocn_vars": ocn_vars,
        "var_aliases": var_aliases,
        "normalization": "palm_gda_zscore",
        "normalization_stats_file": out_stats.name,
        "normalization_stats_path": str(out_stats),
        "handle_nan": True,
        "nan_fill_value": 0.0,
        "atm_shape": atm_shape,
        "ocn_shape": ocn_shape,
        "vert_mode": "fold_z",
        "streaming": True,
        "parallel_backend": "process" if preprocess_workers > 1 and parallel_backend in {"auto", "process", "multiprocessing", "local"} else "single",
        "preprocess_workers": preprocess_workers,
        "time_chunk": time_chunk,
        "num_source_pairs": len(plan),
        "num_frames": len(frame_files),
        "frame_files": frame_files,
        "time_values": list(time_values),
        "aggregate_file": None,
    }
    with open(out_meta, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    with open(output_dir / "_frame_manifest.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    return None, out_meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_DIR / "configs" / "preprocess_nc.json"))
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--parallel-backend", choices=["auto", "mpi", "mpi4py", "process", "multiprocessing", "local", "single"], default=None)
    parser.add_argument("--preprocess-workers", type=int, default=None)
    args = parser.parse_args()
    with open(args.config, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    if args.parallel_backend is not None:
        cfg["parallel_backend"] = args.parallel_backend
    if args.preprocess_workers is not None:
        cfg["preprocess_workers"] = args.preprocess_workers
    out_pt, out_meta = preprocess(cfg, args.output_dir or cfg.get("output_dir", DATA_DIR))
    if out_pt is not None:
        print(f"saved {out_pt}")
    print(f"saved {out_meta}")


if __name__ == "__main__":
    main()
