import re
from tqdm import tqdm



class DataFrameProgress(object):
    def __init__(self, **kwds):
        self._bar = tqdm(**kwds)
        self._data = None

    def __call__(self, state, df):
        self.update_data(df)
        self.update_progress_bar(len(df))
        return self, df

    def update_data(self, df):
        raise NotImplementedError

    def update_progress_bar(self, num_rows):
        self._bar.update(num_rows)
        self._bar.set_postfix(self._data.to_dict())

    def close(self):
        #self._bar.flush()
        self._bar.close()

    def save(self, path):
        self._data.to_csv(path, index=False)



def kmg_bases_to_int(value: str) -> int:
    try:
        result = int(value)
    except:
        result = None
    if result is not None:
        return result

    value_re = re.compile(r"(\d+)([KkMmGg])b*")
    m = value_re.match(value)
    if not m:
        raise ValueError(f"Invalid string: {value}")
    value = int(m.group(1))
    exponent = {"k": 1e3, "m": 1e6, "g": 1e9}[m.group(2).lower()]
    return value * int(exponent)
