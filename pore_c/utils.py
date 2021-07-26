import re
from logging import getLogger

import numpy as np
from tqdm import tqdm

logger = getLogger(__name__)


PHRED_TO_PROB = np.power(10, (np.arange(256, dtype=float) / -10.0))


def mean_qscore(quals):
    return -10 * np.log10(PHRED_TO_PROB[quals].mean())


class DataFrameProgress:
    def __init__(self, save_to=None, **kwds):
        self._bar = tqdm(**kwds)
        self._data = None
        self._save_to = save_to

    def __call__(self, _, df):
        self.update_data(df)
        self.update_progress_bar(len(df))
        return self, df

    def update_data(self, df):
        raise NotImplementedError

    def update_progress_bar(self, num_rows):
        self._bar.update(num_rows)
        self._bar.set_postfix(self._data.to_dict())

    def close(self):
        # self._bar.flush()
        self._bar.close()
        if self._save_to:
            self.save(self._save_to)

    def save(self, path):
        self._data.to_csv(path, index=False)


def kmg_bases_to_int(value: str) -> int:
    try:
        result = int(value)
    except Exception as _:  # noqa
        result = None
    if result is not None:
        return result

    value_re = re.compile(r"(\d+(?:\.\d+)?)([KkMmGg])([bB])?")
    m = value_re.match(value.strip())
    if not m:
        raise ValueError(f"Invalid string: {value}")

    value = float(m.group(1))
    exponent = {"k": 1e3, "m": 1e6, "g": 1e9}[m.group(2).lower()]
    return value * int(exponent)
