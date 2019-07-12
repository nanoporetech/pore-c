from itertools import combinations
from typing import Dict, Union
from pore_c.model import AlignDf, Chrom, PairDf
from pore_c.io import PairFileWriter
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys



def convert_align_df_to_pairs(align_df: AlignDf, chrom_lengths: Dict[Chrom, int], genome_assembly: str, pair_file: str, n_workers: int = 1):
    from streamz import Stream
    from time import sleep
    from dask.distributed import Client, LocalCluster

    parallel = n_workers > 1
    if parallel:
        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)
    writer = PairFileWriter(pair_file, chrom_lengths, genome_assembly, columns = list(PairDf.DTYPE.keys()))


    def pairs_progress(state, pair_df, progress_bar=None):
        pair_counts = pair_df.pair_type.value_counts().to_dict()
        pair_counts['total'] = len(pair_df)
        if state is None:
            state = pair_counts
        else:
            # print(state, batch_stats)
            for key, val in pair_counts.items():
                state[key] += val
        if progress_bar is not None:
            progress_bar.update(len(pair_df))
            progress_bar.set_postfix(state)
        return state, pair_df


    batch_progress_bar = tqdm(total=align_df.npartitions, desc="Batches submitted: ", unit=" batches", position=0)
    pair_progress_bar = tqdm(total=None, desc="Pairs written: ", unit=" pairs", position=1)

    df_stream = Stream()
    # stream that holds the filtered/processed alignments
    if parallel:
        pair_stream = (
            df_stream
            .scatter()
            .map(to_pairs)
            .buffer(n_workers)
            .gather()
        )
    else:
        pair_stream = df_stream.map(to_pairs)

    write_sink = (
        pair_stream
        .accumulate(pairs_progress, returns_state=True, start=None, progress_bar=pair_progress_bar)
        .sink(writer)
    )

    use_cols = ['pass_filter', 'read_idx', 'read_name', 'read_start', 'read_end', 'chrom', 'start',  'strand', 'fragment_id', 'align_idx', 'fragment_start', 'fragment_end']
    for partition in range(align_df.npartitions):
        _df = align_df.partitions[partition].loc[:, use_cols]
        df_stream.emit(_df.compute()) #.head(1000))
        batch_progress_bar.update(1)

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                break
        client.close()
        cluster.close()

    batch_progress_bar.close()
    pair_progress_bar.close()
    sys.stdout.write('\n\n')


def to_pairs(df):
    res = []
    keep_segments = (
        df.query("pass_filter == True")
        .replace({'strand': {True: '+', False: '-'}})
        .astype({'strand': PairDf.DTYPE['strand1']})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue
        read_id = "{}:{}".format(read_df.read_name.values[0], read_idx)
        rows = list(
            read_df
            .sort_values(['read_start'], ascending=True)
            .assign(
                pos_on_read = lambda x: np.arange(len(x)),
                fragment_midpoint = lambda x: np.rint((x.fragment_start + x.fragment_end) * 0.5).astype(int)
            )
            .itertuples()
        )
        for (pos1, pos2) in combinations(rows, 2):
            if pos1.align_idx == pos2.align_idx:
                raise ValueError
            res.append(position_pair_to_tuple(pos1, pos2, read_id))
    return pd.DataFrame(res, columns=PairDf.DTYPE.keys()).astype(PairDf.DTYPE)


def position_pair_to_tuple(pos1, pos2, read_id):
    switch_order = False #reorder to make uppper triangle
    distance_on_read = pos2.read_start - pos1.read_end
    direct_neighbors = (pos2.pos_on_read - pos1.pos_on_read) == 1
    if pos1.chrom == pos2.chrom:
        if pos1.fragment_midpoint > pos2.fragment_midpoint:
            switch_order = True
    elif pos1.chrom > pos2.chrom:
        switch_order = True
    if switch_order:
        pos1, pos2 = pos2, pos1
    return (
         read_id,
         pos1.chrom,
         pos1.fragment_midpoint,
         pos2.chrom,
         pos2.fragment_midpoint,
         pos1.strand,
         pos2.strand,
         'DJ' if direct_neighbors else 'IJ',
         pos1.fragment_id,
         pos2.fragment_id,
         pos1.align_idx,
         pos2.align_idx,
         distance_on_read
    )
