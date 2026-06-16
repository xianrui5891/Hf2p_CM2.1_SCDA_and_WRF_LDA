import numpy as np
import os
import json
import time

NMC_BACKGROUND_COVARIANCE_PATH = '/data/ouc/zhuxianrui/vae-pywrf/nmc_test/nmc_outputs_q11_12/latent_nmc_background_covariance.npz'
TIMING_LOG_PATH = os.environ.get("WRF_LDA_TIMING_LOG_PATH", "./wrf_lda_timing.jsonl")

class transpond_class:
    def __init__(self):
        self.set_up_flg = False
        self.save_flg = True
        self.comm = None
        self.rank = 0
        self.size = 1
        self.timing_calls = 0
        self.timing_keys = [
            "callback_total",
            "da_compute",
        ]
        self.timing_totals_max = {key: 0.0 for key in self.timing_keys}
        self.timing_totals_avg = {key: 0.0 for key in self.timing_keys}

    def _set_up(self, lon, lat):
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        except :
            # If no mpi4py, assume single process
            self.comm = None
            self.rank = 0
            self.size = 1

        if self.rank == 0:
            import torch
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
            from vaelda.vae_lda import latent_space_da
            self.lda_model = latent_space_da(
                dev=self.device,
                obs_path='./obs',
                vae_model_path='./vaelda/best_model.pth',
                nmc_background_covariance_path=NMC_BACKGROUND_COVARIANCE_PATH,
            )
        else:
            self.lda_model = None
        
        self.set_up_flg = True

        if self.rank == 0:
            from utils.data_output_utils import NCWriter, FieldPlotter
            self.saver = NCWriter("./python_nc_data", lon, lat)
            self.plotter = FieldPlotter("./python_plots", lon, lat)
            from utils.obsdata_utils import obs_data_loader
            self.obs_data = obs_data_loader(obs_dir='./obs').load_all()
        else:
            self.saver = None
            self.plotter = None
            self.obs_data = None

    def _rank0_print(self, message: str):
        if self.rank == 0:
            print(message)

    def _write_timing_jsonl(self, payload: dict):
        if self.rank != 0:
            return
        log_dir = os.path.dirname(TIMING_LOG_PATH)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        with open(TIMING_LOG_PATH, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, ensure_ascii=False) + "\n")

    def _reduce_timing(self, timings: dict) -> dict:
        reduced = {}
        for key in self.timing_keys:
            value = float(timings.get(key, 0.0))
            if self.comm is not None:
                max_value = self.comm.allreduce(value, op=self._mpi_max())
                sum_value = self.comm.allreduce(value, op=self._mpi_sum())
                avg_value = sum_value / max(float(self.size), 1.0)
            else:
                max_value = value
                avg_value = value
            reduced[key] = {"max": max_value, "avg": avg_value}
        return reduced

    def _mpi_max(self):
        from mpi4py import MPI
        return MPI.MAX

    def _mpi_sum(self):
        from mpi4py import MPI
        return MPI.SUM

    def _record_timing(self, time_str: str, timings: dict, shapes: dict):
        reduced = self._reduce_timing(timings)
        if self.rank != 0:
            return

        self.timing_calls += 1
        for key in self.timing_keys:
            self.timing_totals_max[key] += reduced[key]["max"]
            self.timing_totals_avg[key] += reduced[key]["avg"]

        payload = {
            "event": "wrf_lda_callback_timing",
            "time_str": time_str,
            "rank_count": self.size,
            "shapes": shapes,
            "timings": reduced,
            "note": "max is the parallel wall-clock cost; avg is mean rank-local cost.",
        }
        self._write_timing_jsonl(payload)

        print(
            "[WRF_LDA_PY_TIMING] "
            f"time={time_str} ranks={self.size} "
            f"callback max/avg={reduced['callback_total']['max']:.4f}/{reduced['callback_total']['avg']:.4f}s "
            f"da max/avg={reduced['da_compute']['max']:.4f}/{reduced['da_compute']['avg']:.4f}s"
        )

    def print_timing_summary(self):
        if self.rank != 0 or self.timing_calls == 0:
            return
        print("=== WRF LDA Python timing summary ===")
        print(f"calls: {self.timing_calls}, ranks: {self.size}")
        for key in self.timing_keys:
            print(
                f"{key}: max_sum={self.timing_totals_max[key]:.4f}s, "
                f"avg_rank_sum={self.timing_totals_avg[key]:.4f}s"
            )
        print(f"jsonl: {TIMING_LOG_PATH}")

    def _save_data_(self, time_str, u_2d, v_2d, t_2d, q_2d, pts, data_type):
        self.saver.save_2d(time_str=time_str, var_name='u', data2d=u_2d, data_type=data_type, units='m/s', long_name='zonal wind at level 0')
        self.saver.save_2d(time_str=time_str, var_name='v', data2d=v_2d, data_type=data_type, units='m/s', long_name='meridional wind at level 0')
        self.saver.save_2d(time_str=time_str, var_name='t', data2d=t_2d, data_type=data_type, units='K', long_name='temperature at level 0')
        self.saver.save_2d(time_str=time_str, var_name='q', data2d=q_2d, data_type=data_type, units='kg/kg', long_name='humidity')
        self.plotter.plot_and_save(time_str=time_str, var_name='u', data2d=u_2d, data_type=data_type, points=pts)
        self.plotter.plot_and_save(time_str=time_str, var_name='v', data2d=v_2d, data_type=data_type, points=pts)
        self.plotter.plot_and_save(time_str=time_str, var_name='t', data2d=t_2d, data_type=data_type, points=pts)
        self.plotter.plot_and_save(time_str=time_str, var_name='q', data2d=q_2d, data_type=data_type, points=pts)

    def __call__(self,
            time_str : str,
            u : np.array,
            v : np.array,
            t : np.array,
            q : np.array,
            xlat: np.array,
            xlon: np.array,
            x_size : int = None,
            y_size : int = None,
            z_size : int = None
    ):
        """
            This function should be invoked by the fortran via f2py callback.
            The data from WRF will be transmitted into this function,
            and then go to da process.
            -------------------------------------------
            time_str: the time string (current_timestr from fortran)
            u/v/t/q: rank-2 arrays with bounds (x_size, z_size) - MODIFY IN PLACE

            Note: Arrays are passed by reference - modify directly, don't return copies!
            Xianrui Zhu 2025/6/8
        """
        
        # Get MPI rank for log optimization
        callback_t0 = time.perf_counter()
        timings = {key: 0.0 for key in self.timing_keys}
        time_str_copy = time_str.strip()
        if not self.set_up_flg:
            self._set_up(xlon.copy(), xlat.copy())

        print(f"python got the data of {time_str_copy}")

        # Check if all arrays have consistent shapes
        
        #print(f"u shape: {u.shape}, range: [{u.min():.6f}, {u.max():.6f}]")
        #print(f"v shape: {v.shape}, range: [{v.min():.6f}, {v.max():.6f}]")
        #print(f"t shape: {t.shape}, range: [{t.min():.6f}, {t.max():.6f}]")
        #print(f"q shape: {q.shape}, range: [{q.min():.6f}, {q.max():.6f}]")
        #print(f"xlon shape: {xlat.shape}, range: [{xlat.min():.6f}, {xlat.max():.6f}]")
        #print(f"xlat shape: {xlon.shape}, range: [{xlon.min():.6f}, {xlon.max():.6f}]")
        
        # Only print full array contents on rank 0 to avoid large log files
        
        obs_coords = []
        if self.rank == 0:
            #print("u: {}\nv: {}\nt: {}\nq:{}\nlat:{}\nlon:{}".format(u, v, t, q, xlat[0,:], xlon[:,0]))

            if time_str_copy in self.obs_data:
                for coor, dicts in self.obs_data[time_str_copy].items():
                    y, x = coor
                    if xlon[:,0][0] < x < xlon[:,0][-1] and xlat[0,:][0] < y < xlat[0,:][-1]:
                        obs_coords.append((x, y))
                        
            if(self.save_flg):
                self._save_data_(time_str_copy, u[:,0,:].copy(), v[:,0,:].copy(), t[:,0,:].copy(), q.copy(), obs_coords, "original")
        else:
            print("Array contents logged only on rank 0 to save space")
        
        # Perform data assimilation processing here
        # Note: modify u, v, t directly, don't assign to new variables
        # Example: u[:] = some_processing(u)

        if self.rank == 0:
            input_fields = [u[:,0,:].copy(), v[:,0,:].copy(), t[:,0,:].copy(), q.copy()]
            axis = [xlon[:,0].copy(), xlat[0,:].copy()]

            da_t0 = time.perf_counter()
            decoded = self.lda_model.do_da(time_str_copy, input_fields, axis) # 先x后y
            timings["da_compute"] = time.perf_counter() - da_t0

            du, dv, dt, dq = decoded
            du = du - u[:,0,:]
            dv = dv - v[:,0,:]
            dt = dt - t[:,0,:]
            dq = dq - q
            u[:,0,:] = u[:,0,:] + 0.3*du
            v[:,0,:] = v[:,0,:] + 0.3*dv
            t[:,0,:] = t[:,0,:] + 0.3*dt
            if(self.save_flg):
                self._save_data_(time_str_copy, u[:,0,:].copy(), v[:,0,:].copy(), t[:,0,:].copy(), q.copy(), obs_coords, "da")
                self._save_data_(time_str_copy, du, dv, dt, dq, obs_coords, "increment")

            #print(f"increment:\ndelta_u:{du}\ndelta_v:{dv}\ndelta_t:{dt}\ndelta_q:{dq}")
        
        if self.comm is not None:
            self.comm.Bcast(u)
            self.comm.Bcast(v)
            self.comm.Bcast(t)
        timings["callback_total"] = time.perf_counter() - callback_t0

        shapes = {
            "u": tuple(int(v) for v in u.shape),
            "v": tuple(int(v) for v in v.shape),
            "t": tuple(int(v) for v in t.shape),
            "q": tuple(int(v) for v in q.shape),
            "xlat": tuple(int(v) for v in xlat.shape),
            "xlon": tuple(int(v) for v in xlon.shape),
        }
        self._record_timing(time_str_copy, timings, shapes)
        print(f"Data assimilation processing completed for {time_str_copy}")
            
        return
