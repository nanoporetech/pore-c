import pytest
import pandas as pd
from pore_c.model import GenomeIntervalDf
pd.options.display.width = 1200
pd.options.display.expand_frame_repr = True
pd.options.display.max_colwidth = 200
import numpy as np


def _to_interval_df(d):
    return pd.DataFrame(d, columns=['chrom', 'start', 'end', 'id'])

@pytest.mark.parametrize(
    "query, target, overlaps",
    [
        (
            [('chr1', 0, 10, 0),], # query in form chrom, start, end, id
            [('chr1', 0, 10, 0),], # target in form chrom, start, end, id
            [(0, 0, 10)] # overlap in form query_idx, target_idx, bases overlap
        ),
        (
            [('chr1', 0, 10, 0),],
            [('chr2', 0, 10, 0),],
            []
        ),
        (
            [('chr1', 0, 10, 0),],
            [('chr1', 10, 11, 0),],
            []
        ),
        (
            [('chr1', 0, 10, 0),],
            [('chr1', 9, 10, 0),],
            [(0, 0, 1)]
        ),
        (
            [('chr1', 0, 10, 0),],
            [('chr1', 0, 1, 0),],
            [(0, 0, 1)]
        ),
        (
            [('chr1', 0, 10, 0),],
            [('chr1', 0, 1, 0), ('chr1', 5, 10, 1)],
            [(0, 0, 1), (0, 1, 5)]
        ),
        (
            [('chr1', 0, 1, 0), ('chr1', 5, 10, 1)],
            [('chr1', 0, 10, 0),],
            [(0, 0, 1), (1, 0, 5)]
        ),
    ]
)
def test_overlaps(query, target, overlaps):
    q = _to_interval_df(query)
    t = _to_interval_df(target)
    olap_df = q.ginterval.overlap(t)
    if olap_df is None:
        res = set([])
    else:
        res = set((tuple(v) for v in olap_df[['index', 'other_index', 'overlap_length']].values))
    assert(res == set(overlaps))

