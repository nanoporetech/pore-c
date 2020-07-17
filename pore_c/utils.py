import re
from contextlib import AbstractContextManager
from logging import getLogger
from time import sleep
from typing import Optional

import numpy as np
from dask.distributed import Client, LocalCluster
from tqdm import tqdm


logger = getLogger(__name__)


PHRED_TO_PROB = np.power(10, (np.arange(256, dtype=float) / -10.0))


def mean_qscore(quals):
    return -10 * np.log10(PHRED_TO_PROB[quals].mean())


class DaskExecEnv(AbstractContextManager):
    def __init__(
        self,
        n_workers: int = 1,
        processes: bool = True,
        threads_per_worker: int = 1,
        scheduler_port: int = 0,
        dashboard_port: Optional[int] = None,
    ):
        self._cluster_kwds = {
            "processes": processes,
            "n_workers": n_workers,
            "scheduler_port": scheduler_port,
            "dashboard_address": f"127.0.0.1:{dashboard_port}",
            "threads_per_worker": threads_per_worker,
        }
        if dashboard_port is None:
            self._cluster_kwds["dashboard_address"] = None
        self._cluster, self._client = None, None

    def scatter(self, data):
        return self._client.scatter(data)

    def __enter__(self):
        self._cluster = LocalCluster(**self._cluster_kwds)
        self._client = Client(self._cluster)
        logger.debug(f"Cluster started: {self._cluster}")
        logger.debug(f"Client started: {self._client}")
        return self

    def __exit__(self, *args):
        if self._cluster:
            max_tries = 10
            backoff = 2
            delay = 1
            while max_tries > 1:
                processing = self._client.processing()
                still_running = [len(v) > 0 for k, v in processing.items()]
                if any(still_running):
                    sleep(delay)
                    max_tries -= 1
                    delay = delay * backoff
                else:
                    sleep(1)
                    break
            self._client.close()
            self._cluster.close()


class DataFrameProgress(object):
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

    value_re = re.compile(r"(\d+)([KkMmGg])b*")
    m = value_re.match(value)
    if not m:
        raise ValueError(f"Invalid string: {value}")
    value = int(m.group(1))
    exponent = {"k": 1e3, "m": 1e6, "g": 1e9}[m.group(2).lower()]
    return value * int(exponent)
