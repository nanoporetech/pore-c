import re
from contextlib import AbstractContextManager
from time import sleep

from dask.distributed import Client, LocalCluster
from tqdm import tqdm


class DaskExecEnv(AbstractContextManager):
    def __init__(self, n_workers: int = 1, empty_queue=False):
        self.parallel = n_workers > 1
        self.n_workers = n_workers
        self.empty_queue = empty_queue

    def scatter(self, data):
        return self._client.scatter(data)

    def __enter__(self):
        if self.parallel:
            self._cluster = LocalCluster(processes=True, n_workers=self.n_workers, threads_per_worker=1)
            self._client = Client(self._cluster)
        else:
            self._cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=1)
            self._client = Client(self._cluster)
        return self

    def __exit__(self, *args):
        if self.parallel and self.empty_queue:
            while True:
                processing = self._client.processing()
                still_running = [len(v) > 0 for k, v in processing.items()]
                if any(still_running):
                    sleep(5)
                else:
                    sleep(1)
                    break
            self._client.close()
            self._cluster.close()


class DataFrameProgress(object):
    def __init__(self, save_to=None, **kwds):
        self._bar = tqdm(**kwds)
        self._data = None
        self._save_to = None

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

    value_re = re.compile(r"(\d+)([KkMmGg])b*")
    m = value_re.match(value)
    if not m:
        raise ValueError(f"Invalid string: {value}")
    value = int(m.group(1))
    exponent = {"k": 1e3, "m": 1e6, "g": 1e9}[m.group(2).lower()]
    return value * int(exponent)
