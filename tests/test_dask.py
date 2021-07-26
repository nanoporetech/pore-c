import dask.dataframe as dd
import pandas as pd
import pytest

from pore_c.dask import ExecContext, SchedulerType


@pytest.fixture(
    params=[
        (SchedulerType.default, {}),
        (SchedulerType.synchronous, {}),
        (SchedulerType.threads, {"n_workers": 10}),
        (SchedulerType.processes, {"n_workers": 10}),
        (SchedulerType.local_cluster, {"processes": False, "n_workers": 10, "threads_per_worker": 1}),
    ],
    ids=["sched-default", "sched-sync", "sched-threads", "sched-processes", "sched-local_cluster"],
)
def exec_ctx(request):
    return ExecContext(request.param[0], request.param[1])


def test_dask_env(exec_ctx):

    with exec_ctx:
        df = dd.from_pandas(pd.DataFrame({"a": range(100), "b": range(100)}), npartitions=10)
        assert df.sum().sum().compute() == 9900
