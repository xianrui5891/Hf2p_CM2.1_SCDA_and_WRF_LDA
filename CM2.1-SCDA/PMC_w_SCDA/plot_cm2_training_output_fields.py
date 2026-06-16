#!/usr/bin/env python3
"""
Plot CM2 WCDA/SCDA fields and their differences.

The default paths follow plot_obs_model_error.py:
  wcda: /data/assim-f2py-test/work_ctl/cm2_wcda/{atmos,ocean}_3hourly.nc
  scda: /data/assim-f2py-test/work_ctl/cm2_scda_final7/{atmos,ocean}_3hourly.nc

Each variable/time/case is saved as a separate image. Color limits are scanned
once over the selected real-time interval, then reused for every image of that
variable. Difference images use separately scanned fixed limits.
"""

from __future__ import annotations

import argparse
import ast
import json
import logging
import sys
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import xarray as xr


PROJECT_DIR = Path(__file__).resolve().parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from data_preprocess.preprocess_nc import _prepare_var
from plot_obs_model_error import FILES, WORK_DIRS_DEFAULT
from utils.plotting import FieldPlotter


DEFAULT_PREPROCESS_CONFIG = PROJECT_DIR / "configs" / "preprocess_nc.json"
DEFAULT_OUTDIR = "./plots/cm2_wcda_scda_fields_19821026_19821113"
DEFAULT_START = "1982-10-26"
DEFAULT_END = "1982-11-13"
DEFAULT_CASES = "wcda,scda"

VARIABLES = {
    "u": {"domain": "atm", "channel": "UCOMP", "label": "u"},
    "v": {"domain": "atm", "channel": "VCOMP", "label": "v"},
    "t": {"domain": "atm", "channel": "TEMP", "label": "t"},
    "p": {"domain": "atm", "channel": "PS", "label": "p"},
    "ssu": {"domain": "ocn", "channel": "U_SURF", "label": "ssu"},
    "ssv": {"domain": "ocn", "channel": "V_SURF", "label": "ssv"},
    "sst": {"domain": "ocn", "channel": "SST", "label": "sst"},
    "ssh": {"domain": "ocn", "channel": "ETA_T", "label": "ssh"},
}


class LimitStats:
    def __init__(self) -> None:
        self.lo = np.inf
        self.hi = -np.inf

    def update(self, values: np.ndarray) -> None:
        if np.ma.isMaskedArray(values):
            finite = np.asarray(values.compressed(), dtype=np.float64)
        else:
            finite = np.asarray(values, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        self.lo = min(self.lo, float(finite.min()))
        self.hi = max(self.hi, float(finite.max()))

    def as_limits(self) -> tuple[float | None, float | None]:
        if not np.isfinite(self.lo) or not np.isfinite(self.hi):
            return None, None
        if self.lo == self.hi:
            delta = max(abs(self.lo) * 1e-6, 1e-12)
            return self.lo - delta, self.hi + delta
        return self.lo, self.hi

    def as_symmetric_limits(self) -> tuple[float | None, float | None]:
        if not np.isfinite(self.lo) or not np.isfinite(self.hi):
            return None, None
        bound = max(abs(float(self.lo)), abs(float(self.hi)), 1e-12)
        return -bound, bound


def setup_logging(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    stream = logging.StreamHandler()
    stream.setFormatter(fmt)
    file_handler = logging.FileHandler(outdir / "plot_cm2_training_output_fields.log")
    file_handler.setFormatter(fmt)
    logger.addHandler(stream)
    logger.addHandler(file_handler)


def parse_time_window(start: str, end: str) -> tuple[pd.Timestamp, pd.Timestamp]:
    start_ts = pd.to_datetime(start)
    end_ts = pd.to_datetime(end)
    if end_ts.time() == pd.Timestamp(end_ts.date()).time():
        end_ts = end_ts + pd.Timedelta(days=1) - pd.Timedelta(nanoseconds=1)
    if start_ts > end_ts:
        raise ValueError(f"start time {start_ts} is after end time {end_ts}")
    return start_ts, end_ts


def normalize_time(value) -> pd.Timestamp:
    if isinstance(value, pd.Timestamp):
        return value
    try:
        return pd.to_datetime(value)
    except Exception:
        pass
    try:
        return pd.Timestamp(
            year=int(getattr(value, "year")),
            month=int(getattr(value, "month", 1)),
            day=int(getattr(value, "day", 1)),
            hour=int(getattr(value, "hour", 0)),
            minute=int(getattr(value, "minute", 0)),
            second=int(getattr(value, "second", 0)),
        )
    except Exception:
        return pd.to_datetime(str(value))


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_preprocess_config(path: str | Path | None) -> dict:
    if not path:
        return {}
    config_path = Path(path)
    if not config_path.exists():
        logging.warning("Preprocess config not found: %s; using built-in defaults.", config_path)
        return {}
    return load_json(config_path)


def open_raw_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=True, decode_timedelta=True, mask_and_scale=True, cache=False)


def raw_file_for_domain(workdir: Path, domain: str) -> Path:
    return workdir / FILES[domain]


def cm2_output_paths(workdir: Path) -> dict[str, list[Path]]:
    paths = {domain: [raw_file_for_domain(workdir, domain)] for domain in ("atm", "ocn")}
    for domain, domain_paths in paths.items():
        for path in domain_paths:
            if not path.exists():
                raise FileNotFoundError(f"Missing CM2 {domain} output file: {path}")
    return paths


def parse_cases(value: str) -> list[str]:
    cases = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not cases:
        raise ValueError("At least one case is required.")
    unknown = [case for case in cases if case not in WORK_DIRS_DEFAULT]
    if unknown:
        raise ValueError(f"Unknown cases: {unknown}; choices are {sorted(WORK_DIRS_DEFAULT)}")
    return cases


def parse_workdirs(value: str | dict | None) -> dict[str, str]:
    if value is None:
        return dict(WORK_DIRS_DEFAULT)
    if isinstance(value, dict):
        return dict(value)
    parsed = ast.literal_eval(value)
    if not isinstance(parsed, dict):
        raise ValueError("--workdirs must be a Python dict literal.")
    return {str(key): str(path) for key, path in parsed.items()}


def nc_field_arrays(
    paths_by_domain: dict[str, list[Path]],
    start: pd.Timestamp,
    end: pd.Timestamp,
    display_vars: Iterable[str],
    *,
    var_aliases: dict,
    vert_select: dict,
) -> Iterable[tuple[str, pd.Timestamp, np.ndarray]]:
    for domain in ("atm", "ocn"):
        requested = [name for name in display_vars if VARIABLES[name]["domain"] == domain]
        if not requested:
            continue
        for path in paths_by_domain.get(domain, []):
            ds = open_raw_dataset(path)
            try:
                prepared: dict[str, xr.DataArray] = {}
                for display_var in requested:
                    channel = VARIABLES[display_var]["channel"]
                    prepared[display_var] = _prepare_var(
                        ds,
                        channel,
                        domain,
                        vert_select=vert_select,
                        var_aliases=var_aliases,
                    )

                for display_var, da in prepared.items():
                    time_dim = da.dims[0]
                    raw_times = [normalize_time(value) for value in np.asarray(da[time_dim].values)]
                    selected = [i for i, timestamp in enumerate(raw_times) if start <= timestamp <= end]
                    for index in selected:
                        timestamp = raw_times[index]
                        data = np.asarray(da.isel({time_dim: index}).values, dtype=np.float32)
                        yield display_var, timestamp, np.ma.masked_invalid(data)
            finally:
                ds.close()


def time_index_map(da: xr.DataArray) -> dict[pd.Timestamp, int]:
    time_dim = da.dims[0]
    out = {}
    for index, value in enumerate(np.asarray(da[time_dim].values)):
        out[normalize_time(value)] = index
    return out


def diff_field_arrays(
    scda_paths: dict[str, list[Path]],
    wcda_paths: dict[str, list[Path]],
    start: pd.Timestamp,
    end: pd.Timestamp,
    display_vars: Iterable[str],
    *,
    var_aliases: dict,
    vert_select: dict,
) -> Iterable[tuple[str, pd.Timestamp, np.ndarray]]:
    for domain in ("atm", "ocn"):
        requested = [name for name in display_vars if VARIABLES[name]["domain"] == domain]
        if not requested:
            continue
        scda_domain_paths = scda_paths.get(domain, [])
        wcda_domain_paths = wcda_paths.get(domain, [])
        if len(scda_domain_paths) != 1 or len(wcda_domain_paths) != 1:
            raise ValueError("Difference output expects exactly one NetCDF file per case/domain.")

        ds_scda = open_raw_dataset(scda_domain_paths[0])
        ds_wcda = open_raw_dataset(wcda_domain_paths[0])
        try:
            for display_var in requested:
                channel = VARIABLES[display_var]["channel"]
                da_scda = _prepare_var(
                    ds_scda,
                    channel,
                    domain,
                    vert_select=vert_select,
                    var_aliases=var_aliases,
                )
                da_wcda = _prepare_var(
                    ds_wcda,
                    channel,
                    domain,
                    vert_select=vert_select,
                    var_aliases=var_aliases,
                )
                scda_indices = time_index_map(da_scda)
                wcda_indices = time_index_map(da_wcda)
                common_times = sorted(timestamp for timestamp in scda_indices if timestamp in wcda_indices and start <= timestamp <= end)
                for timestamp in common_times:
                    scda_data = np.asarray(da_scda.isel({da_scda.dims[0]: scda_indices[timestamp]}).values, dtype=np.float32)
                    wcda_data = np.asarray(da_wcda.isel({da_wcda.dims[0]: wcda_indices[timestamp]}).values, dtype=np.float32)
                    yield display_var, timestamp, np.ma.masked_invalid(scda_data - wcda_data)
        finally:
            ds_scda.close()
            ds_wcda.close()


def update_limits(
    limits: dict[str, LimitStats],
    arrays: Iterable[tuple[str, pd.Timestamp, np.ndarray]],
    source_name: str,
) -> int:
    count = 0
    for display_var, _, arr in arrays:
        limits[display_var].update(arr)
        count += 1
    logging.info("Scanned %d fields for fixed color limits from %s", count, source_name)
    return count


def time_token(timestamp: pd.Timestamp) -> str:
    return pd.to_datetime(timestamp).strftime("%Y%m%dT%H%M%S")


def orient_for_field_plotter(arr: np.ndarray) -> np.ndarray:
    """FieldPlotter transposes input internally; pre-transpose NetCDF y/x data."""
    masked = np.ma.masked_invalid(arr)
    if masked.ndim != 2:
        raise ValueError("data2d must be 2D")
    return masked.T


def plot_arrays(
    outdir: Path,
    source_name: str,
    limits: dict[str, LimitStats],
    arrays: Iterable[tuple[str, pd.Timestamp, np.ndarray]],
    *,
    symmetric: bool = False,
) -> int:
    plotter = FieldPlotter(outdir, cmap="bwr", use_runtime_cmap=False)
    count = 0
    for display_var, timestamp, arr in arrays:
        label = VARIABLES[display_var]["label"]
        vmin, vmax = limits[display_var].as_symmetric_limits() if symmetric else limits[display_var].as_limits()
        plotter.plot_and_save(
            time_str=time_token(timestamp),
            var_name=label,
            data2d=orient_for_field_plotter(arr),
            data_type=source_name,
            vmin=vmin,
            vmax=vmax,
        )
        count += 1
        if count % 100 == 0:
            logging.info("Saved %d images for %s", count, source_name)
    logging.info("Saved %d images for %s", count, source_name)
    return count


def parse_vars(value: str) -> list[str]:
    requested = [item.strip().lower() for item in value.split(",") if item.strip()]
    unknown = [item for item in requested if item not in VARIABLES]
    if unknown:
        raise ValueError(f"Unknown variables: {unknown}; choices are {sorted(VARIABLES)}")
    return requested


def main(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir).resolve()
    setup_logging(outdir)
    start, end = parse_time_window(args.start, args.end)
    display_vars = parse_vars(args.vars)
    case_limits = {name: LimitStats() for name in display_vars}
    diff_limits = {name: LimitStats() for name in display_vars}

    logging.info("Output directory: %s", outdir)
    logging.info("Real-time window: %s to %s", start, end)
    logging.info("Variables: %s", ", ".join(display_vars))

    preprocess_config = load_preprocess_config(args.preprocess_config)
    var_aliases = preprocess_config.get("var_aliases", {})
    vert_select = preprocess_config.get("vert_select", {})
    work_dirs = parse_workdirs(args.workdirs)
    cases = parse_cases(args.cases)
    logging.info("Cases: %s", ", ".join(cases))

    case_paths: dict[str, dict[str, list[Path]]] = {}
    for case in cases:
        if case not in work_dirs:
            raise KeyError(f"No workdir configured for case {case!r}.")
        workdir = Path(work_dirs[case])
        case_paths[case] = cm2_output_paths(workdir)
        logging.info("Case %s workdir: %s", case, workdir)
        count = update_limits(
            case_limits,
            nc_field_arrays(
                case_paths[case],
                start,
                end,
                display_vars,
                var_aliases=var_aliases,
                vert_select=vert_select,
            ),
            case,
        )
        if count == 0:
            raise RuntimeError(f"No {case} fields found in {workdir} for {start} to {end}.")

    if args.include_diff:
        if "scda" not in case_paths:
            case_paths["scda"] = cm2_output_paths(Path(work_dirs["scda"]))
        if "wcda" not in case_paths:
            case_paths["wcda"] = cm2_output_paths(Path(work_dirs["wcda"]))
        count = update_limits(
            diff_limits,
            diff_field_arrays(
                case_paths["scda"],
                case_paths["wcda"],
                start,
                end,
                display_vars,
                var_aliases=var_aliases,
                vert_select=vert_select,
            ),
            "scda_minus_wcda",
        )
        if count == 0:
            raise RuntimeError(f"No exact common SCDA/WCDA times found for {start} to {end}.")

    for name in display_vars:
        vmin, vmax = case_limits[name].as_limits()
        logging.info("Fixed case color limits for %s: vmin=%s vmax=%s", name, vmin, vmax)
        if args.include_diff:
            diff_vmin, diff_vmax = diff_limits[name].as_symmetric_limits()
            logging.info("Fixed difference color limits for %s: vmin=%s vmax=%s", name, diff_vmin, diff_vmax)

    for case in cases:
        plot_arrays(
            outdir,
            case,
            case_limits,
            nc_field_arrays(
                case_paths[case],
                start,
                end,
                display_vars,
                var_aliases=var_aliases,
                vert_select=vert_select,
            ),
        )

    if args.include_diff:
        plot_arrays(
            outdir,
            "scda_minus_wcda",
            diff_limits,
            diff_field_arrays(
                case_paths["scda"],
                case_paths["wcda"],
                start,
                end,
                display_vars,
                var_aliases=var_aliases,
                vert_select=vert_select,
            ),
            symmetric=True,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frame-dir",
        type=str,
        default=None,
        help="Ignored compatibility option.",
    )
    parser.add_argument(
        "--preprocess-config",
        type=str,
        default=str(DEFAULT_PREPROCESS_CONFIG),
        help="Config used to get variable aliases and vertical selection.",
    )
    parser.add_argument(
        "--workdir",
        type=str,
        default=None,
        help="Ignored compatibility option; use --workdirs/--cases.",
    )
    parser.add_argument(
        "--workdirs",
        type=str,
        default=str({key: WORK_DIRS_DEFAULT[key] for key in ("wcda", "scda")}),
        help="Python dict literal mapping case names to CM2 work directories.",
    )
    parser.add_argument("--cases", type=str, default=DEFAULT_CASES, help="Comma-separated cases to plot, default wcda,scda.")
    parser.add_argument("--outdir", type=str, default=DEFAULT_OUTDIR)
    parser.add_argument("--start", type=str, default=DEFAULT_START, help="Real start time, e.g. 1982-10-26 or 1982-10-26T00:00:00.")
    parser.add_argument("--end", type=str, default=DEFAULT_END, help="Real end time. A date-only value includes the whole day.")
    parser.add_argument("--vars", type=str, default="u,v,t,p,ssu,ssv,sst,ssh")
    parser.add_argument(
        "--sources",
        choices=("training", "cm2_output", "both"),
        default=None,
        help="Ignored compatibility option.",
    )
    parser.add_argument("--include-diff", dest="include_diff", action="store_true", default=True, help="Also plot scda-wcda.")
    parser.add_argument("--no-diff", dest="include_diff", action="store_false", help="Disable scda-wcda plots.")
    main(parser.parse_args())
