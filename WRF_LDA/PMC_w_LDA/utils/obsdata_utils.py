from utils.normalization_utils import var_name, norm_dict
from pathlib import Path
from datetime import datetime

class obs_data_loader:
    def __init__(self, obs_dir='../data/obs',
                  missing=-9999.999,
                  norm_tag=True,
                  keep_fields=None,
                  name_mapping=None):
        self.obs_dir = Path(obs_dir)
        self.missing = float(missing)
        # fields we'll parse (we only keep these ultimately)
        self.keep_fields = keep_fields or ['u_wind','v_wind']
        self.name_mapping = name_mapping or {'u_wind':'U10','v_wind':'V10'}
        self.norm_tag = norm_tag
        self.data = {}
    
    def __is_missing(self, val):
        return abs(float(val) - self.missing) < 1e-6
    
    def __trans_time_fmt(self, raw):
        raw = raw.strip()
        if len(raw) >= 12 and raw[:12].isdigit():
            ts = raw[:12]
            dt = datetime.strptime(ts, '%Y%m%d%H%M')
            return dt.strftime('%Y-%m-%d_%H:%M:00')

    def __read_file__(self, file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [ln.rstrip() for ln in f.readlines() if ln.strip() != '']
        timestamp = self.__trans_time_fmt(lines[0])
        header = lines[1].split()
        header = [h.lower() for h in header]
        idx = {name: i for i, name in enumerate(header)}
        
        required = set([*self.keep_fields,'lat','lon'])
        if not required.issubset(set(header)):
            missing_cols = required - set(header)
            raise KeyError(f'File {file_path} missing required columns: {missing_cols}')

        if timestamp not in self.data:
            self.data[timestamp] = {}

        def normalize(val, val_name):
            standardized = (val - norm_dict[val_name]['mean']) / norm_dict[val_name]['std']
            return standardized

        for raw in lines[2:]:
            parts = raw.split()
            if len(parts) < len(header):
                raise KeyError(f'File {file_path} have misssing data.')

            lat = float(parts[idx['lat']])
            lon = float(parts[idx['lon']])

            key = (lat,lon)
            if key not in self.data[timestamp]:
                self.data[timestamp][key] = {}

            height = float(parts[idx['levm']])
            if height > 100.0:
                continue

            for val_name_obs in self.keep_fields:
                val = float(parts[idx[val_name_obs]])
                val_name = self.name_mapping[val_name_obs]
                if val_name not in self.data[timestamp][key] and not self.__is_missing(float(parts[idx[val_name_obs]])):
                    if self.norm_tag:
                        self.data[timestamp][key][val_name] = normalize(val, val_name)
                    else:
                        self.data[timestamp][key][val_name] = val

        def not_empty(v):
            return v is not None and not (isinstance(v, (str, list, dict, tuple, set)) and len(v) == 0)
        self.data[timestamp] = { k: v for k, v in self.data[timestamp].items() if not_empty(v)}

        return timestamp


    def load_all(self):
        if not self.obs_dir.exists():
            raise FileNotFoundError(f'{self.obs_dir} does not exist.')
        files = sorted([p for p in self.obs_dir.iterdir() if p.is_file()])
        for p in files:
            try:
                self.__read_file__(p)
            except Exception as e:
                # skip problematic files but print a warning
                print(f'Warning: failed to parse {p}: {e}')
        return self.data
    
if __name__ == "__main__":
    test = obs_data_loader(norm_tag=True).load_all()
    print(test)
