import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import ExitStack
from enum import Enum
from typing import Any, Dict

import dask
from dask.distributed import Client, LocalCluster

logger = logging.getLogger(__name__)


class SchedulerType(Enum):
    default = "default"
    synchronous = "synchronous"
    threads = "threads"
    processes = "processes"
    local_cluster = "local_cluster"


class ExecContext(ExitStack):
    def __init__(self, scheduler: SchedulerType, scheduler_kwds: Dict[Any, Any]):
        super().__init__()
        self.scheduler = scheduler
        self.performance_report_path = scheduler_kwds.pop("performance_report", None)
        self.scheduler_kwds = scheduler_kwds

    @classmethod
    def from_cli_opts(cls, **kwds):
        pass

    def __enter__(self):
        super().__enter__()
        logger.debug(f"Setting up env: {self.scheduler} {self.scheduler_kwds}")
        if self.scheduler is SchedulerType.default:
            logger.debug("Scheduler: {self.scheduler}")
        elif self.scheduler is SchedulerType.synchronous:
            self.enter_context(dask.config.set(scheduler="synchronous"))
        elif self.scheduler is SchedulerType.threads:
            n_workers = self.scheduler_kwds.get("n_workers", None)
            if n_workers:
                pool = ThreadPoolExecutor(n_workers)
                self.enter_context(dask.config.set(pool=pool))
            else:
                self.enter_context(dask.config.set(scheduler="threads"))
        elif self.scheduler is SchedulerType.processes:
            n_workers = self.scheduler_kwds.get("n_workers", None)
            if n_workers:
                pool = ProcessPoolExecutor(n_workers)
                self.enter_context(dask.config.set(pool=pool))
            else:
                self.enter_context(dask.config.set(scheduler="processes"))
        elif self.scheduler is SchedulerType.local_cluster:
            cluster = self.enter_context(LocalCluster(**self.scheduler_kwds))
            client = self.enter_context(Client(cluster))
            logger.debug(f"Started LocalCluster:\n{cluster}\n{client}\n")
        return self
