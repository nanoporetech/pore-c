import sys
from itertools import combinations
from logging import getLogger
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from streamz import Stream
from tqdm import tqdm

from pore_c.datasources import IndexedPairFile
from pore_c.io import PairFileWriter
from pore_c.model import (FRAG_IDX_DTYPE, AlignDf, Chrom, GenomeIntervalDf,
                          PairDf, SalsaDf, HicTxtDf)
from pore_c.utils import DataFrameProgress

logger = getLogger(__name__)


class PairsProgress(DataFrameProgress):
    def __init__(self, **kwds):
        kwds["desc"] = "Pairs processed"
        kwds["unit"] = " pairs"
        super(PairsProgress, self).__init__(**kwds)

    def update_data(self, pair_df):
        stats = pair_df["pair_type"].value_counts(dropna=False).to_frame().T.reset_index(drop=True)
        if self._data is None:
            self._data = stats
        else:
            self._data += stats

    def get_summary(self):
        return {key: int(val) for key, val in self._data.iloc[0, :].to_dict().items()}

    def update_progress_bar(self, num_pairs):
        self._bar.update(num_pairs)
        self._bar.set_postfix(self.get_summary())


class MatrixAccumlator(DataFrameProgress):
    def __init__(self, **kwds):
        kwds["desc"] = "Batches processed"
        kwds["unit"] = " batches"
        super(MatrixAccumlator, self).__init__(**kwds)

    def update_data(self, df):
        _df = df.set_index(["bin1_id", "bin2_id"]).sort_index()
        if self._data is None:
            self._data = _df
        else:
            self._data = self._data.add(_df, fill_value=0).astype(int)

    def get_summary(self):
        summary = {
            "num_pixels": len(self._data),
            "max_count": int(self._data["count"].max()),
            "median_count": int(self._data["count"].median()),
            "total_contacts": int(self._data["count"].sum()),
            "diagonal_contacts": int(self._data.query("bin1_id == bin2_id")["count"].sum()),
        }
        return summary

    def update_progress_bar(self, num_bins):
        self._bar.update(1)
        self._bar.set_postfix(self.get_summary())

    def save_coo(self, path):
        if self._data is None:
            raise ValueError(f"No data to write to {path}")
        self._data.to_csv(path, sep="\t", header=None)

    def save_bedgraph(self, path, bin_df):
        raise NotImplementedError


def assign_to_bins(pair_df, bin_df, sort_bins=True):
    _dfs = []
    for pos in ["1", "2"]:
        _df = (
            pair_df.loc[:, [f"chr{pos}", f"pos{pos}"]]
            .rename(columns={f"chr{pos}": "chrom", f"pos{pos}": "start"})
            .assign(end=lambda x: x.start + 1)
            # .assign(start = lambda x: x.start - 1)
        )
        bins = _df.ginterval.assign(bin_df).rename(columns={"other": f"bin{pos}_id"})
        _dfs.append(bins)
    overlaps = pd.concat(_dfs, axis=1)  # .astype(FRAG_IDX_DTYPE)
    has_nulls = overlaps.isna().any(axis=1)
    if has_nulls.any():
        raise ValueError("Some fragments missing overlaps: {}".format(pair_df[has_nulls]))
    overlaps = overlaps.astype(FRAG_IDX_DTYPE)
    if sort_bins:
        # pairs files are in lexographic order, matrxi is in fasta order
        switch_bins = overlaps["bin1_id"] > overlaps["bin2_id"]
        if switch_bins.any():
            logger.warning("Reordering bins")
            overlaps = pd.DataFrame(
                {"bin1_id": overlaps.min(axis=1), "bin2_id": overlaps.max(axis=1)},
                index=overlaps.index,
                dtype=FRAG_IDX_DTYPE,
            )
            switch_bins = overlaps["bin1_id"] > overlaps["bin2_id"]
            assert not switch_bins.any()
    return overlaps


def overlap_count(pair_df, bin_df=None):
    df = (
        assign_to_bins(pair_df, bin_df)
        .groupby(["bin1_id", "bin2_id"])
        .size()
        .rename("count")
        .to_frame()
        .reset_index()
        .astype({"bin1_id": int, "bin2_id": int})
    )
    return df


def convert_pairs_to_matrix(pairs_datasource: IndexedPairFile, resolution: int, coo: Path = None, n_workers: int = 1):
    ds = pairs_datasource
    ds.discover()
    chrom_lengths = ds._chroms
    bin_df = GenomeIntervalDf.fixed_width_bins(chrom_lengths, resolution).set_index("bin_id")
    parallel = n_workers > 1
    if parallel:
        from time import sleep
        from dask.distributed import Client, LocalCluster

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)
        bin_df = client.scatter(bin_df)

    batch_progress_bar = tqdm(total=ds.npartitions, desc="Batches submitted: ", unit=" batches", position=0)
    matrix = MatrixAccumlator(position=1)
    # stream that holds the filtered/processed alignments
    df_stream = Stream()
    if parallel:
        coo_stream = (
            df_stream.scatter()
            .buffer(n_workers)
            .map(overlap_count, bin_df=bin_df)
            .gather()
        )
    else:
        coo_stream = (
            df_stream.map(overlap_count, bin_df=bin_df)
        )

    write_sink = coo_stream.accumulate(matrix, returns_state=True, start=matrix).sink(lambda x: x)  # noqa: F841

    partition_range = range(ds.npartitions)

    for partition in partition_range:
        _df = ds._get_partition(partition, usecols=["chr1", "pos1", "chr2", "pos2"])
        if not _df.empty:
            df_stream.emit(_df)
        else:
            logger.debug(f"Empty partition {partition}")
        batch_progress_bar.update(1)

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                sleep(10)
                break
        client.close()
        cluster.close()

    batch_progress_bar.close()
    matrix.close()
    sys.stderr.write("\n\n")
    if coo:
        logger.info(f"Writing COO to {coo}")
        matrix.save_coo(coo)

    return matrix.get_summary()


def convert_align_df_to_pairs(
    align_df: AlignDf, chrom_lengths: Dict[Chrom, int], genome_assembly: str, pair_file: str, n_workers: int = 1
):
    parallel = n_workers > 1
    if parallel:
        from time import sleep
        from dask.distributed import Client, LocalCluster

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    writer = PairFileWriter(pair_file, chrom_lengths, genome_assembly, columns=list(PairDf.DTYPE.keys()))

    batch_progress_bar = tqdm(total=align_df.npartitions, desc="Batches submitted: ", unit=" batches", position=0)
    pairs_progress = PairsProgress()

    # stream that holds the filtered/processed alignments
    df_stream = Stream()
    if parallel:
        pair_stream = df_stream.scatter().map(to_pairs).buffer(n_workers).gather()
    else:
        pair_stream = df_stream.map(to_pairs)

    write_sink = pair_stream.accumulate(pairs_progress, returns_state=True, start=pairs_progress).sink(  # noqa: F841
        writer
    )

    use_cols = [
        "pass_filter",
        "read_idx",
        "read_name",
        "read_start",
        "read_end",
        "chrom",
        "start",
        "strand",
        "fragment_id",
        "align_idx",
        "fragment_start",
        "fragment_end",
    ]
    for partition in range(align_df.npartitions):
        _df = align_df.partitions[partition].loc[:, use_cols]
        if False:
            df_stream.emit(_df.head(1000))
            if partition > 2:
                break
        else:
            df_stream.emit(_df.compute())
        batch_progress_bar.update(1)

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                sleep(10)
                break
        client.close()
        cluster.close()

    writer.close()
    pairs_progress.close()
    batch_progress_bar.close()
    sys.stderr.write("\n\n")
    return {"contact_types": pairs_progress.get_summary()}


def to_pairs(df: AlignDf) -> PairDf:
    res = []
    keep_segments = (
        df.query("pass_filter == True").replace({"strand": {True: "+", False: "-"}})
        # convert strand to +/- and chrom to string for lexographical sorting to get upper-triangle
        .astype({"strand": PairDf.DTYPE["strand1"], "chrom": str})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue
        read_id = "{}:{}".format(read_df.read_name.values[0], read_idx)
        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .assign(
                pos_on_read=lambda x: np.arange(len(x)),
                fragment_midpoint=lambda x: np.rint((x.fragment_start + x.fragment_end) * 0.5).astype(int),
            )
            .itertuples()
        )
        for (pos1, pos2) in combinations(rows, 2):
            if pos1.align_idx == pos2.align_idx:
                raise ValueError
            res.append(position_pair_to_tuple(pos1, pos2, read_id))
    return pd.DataFrame(res, columns=PairDf.DTYPE.keys()).astype(PairDf.DTYPE)




def position_pair_to_tuple(pos1, pos2, read_id):
    switch_order = False  # reorder to make uppper triangle
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
        "DJ" if direct_neighbors else "IJ",
        pos1.fragment_id,
        pos2.fragment_id,
        pos1.align_idx,
        pos2.align_idx,
        distance_on_read,
    )


def convert_align_df_to_salsa(align_df: AlignDf, salsa_bed: Path, n_workers: int = 1):
    from pore_c.io import SalsaBedFileWriter
    parallel = n_workers > 1
    if parallel:
        from time import sleep
        from dask.distributed import Client, LocalCluster

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    writer = SalsaBedFileWriter(salsa_bed)

    batch_progress_bar = tqdm(total=align_df.npartitions, desc="Batches submitted: ", unit=" batches", position=0)

    # stream that holds the filtered/processed alignments
    df_stream = Stream()
    if parallel:
        pair_stream = df_stream.scatter().map(to_salsa).buffer(n_workers).gather()
    else:
        pair_stream = df_stream.map(to_salsa)

    write_sink = pair_stream.sink(writer)  # noqa: F841
    use_cols = [
        "pass_filter",
        "read_idx",
        "align_idx",
        "read_name",
        "read_start",
        "chrom",
        "start",
        "end",
        "strand",
        "mapping_quality",
    ]
    for partition in range(align_df.npartitions):
        _df = align_df.partitions[partition].loc[:, use_cols]
        df_stream.emit(_df.compute())
        batch_progress_bar.update(1)

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                sleep(10)
                break
        client.close()
        cluster.close()

    num_pairs_written = writer.close()
    batch_progress_bar.close()
    sys.stderr.write("\n\n")
    return num_pairs_written


def to_salsa(df: AlignDf) -> SalsaDf:
    res = []

    keep_segments = (
        df.query("pass_filter == True").replace({"strand": {True: "+", False: "-"}})
        .astype({"strand": SalsaDf.DTYPE['strand'], "chrom": str})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue
        rows = list(read_df.sort_values(["read_start"], ascending=True).itertuples())

        for combi, (pos1, pos2) in enumerate(combinations(rows, 2)):
            if pos1.align_idx == pos2.align_idx:
                raise ValueError
            for rec, paired_id in zip([pos1, pos2], [1, 2]):
                paired_read_id = f"{rec.read_name}_{combi}/{paired_id}"
                res.append(
                    (rec.chrom, rec.start, rec.end , paired_read_id, rec.mapping_quality, rec.strand)
                )
    return pd.DataFrame(res, columns=SalsaDf.DTYPE.keys()).astype(SalsaDf.DTYPE)


def convert_align_df_to_hic(align_df: AlignDf, hic_bed: Path, n_workers: int = 1, max_fragment_id: int = 0):
    from pore_c.io import HicTxtFileWriter
    parallel = n_workers > 1
    if parallel:
        from time import sleep
        from dask.distributed import Client, LocalCluster

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    writer = HicTxtFileWriter(hic_bed)

    batch_progress_bar = tqdm(total=align_df.npartitions, desc="Batches submitted: ", unit=" batch", position=0)

    # stream that holds the filtered/processed alignments
    df_stream = Stream()
    if parallel:
        pair_stream = df_stream.scatter().map(to_hic, max_fragment_id=max_fragment_id).buffer(n_workers).gather()
    else:
        pair_stream = df_stream.map(to_hic, max_fragment_id=max_fragment_id)

    write_sink = pair_stream.sink(writer)  # noqa: F841
    use_cols = [
        "pass_filter",
        "read_idx",
        "align_idx",
        "read_name",
        "read_start",
        "chrom",
        "fragment_id",
        "fragment_end",
        "strand",
        "mapping_quality",
    ]
    for partition in range(align_df.npartitions):
        _df = align_df.partitions[partition].loc[:, use_cols]
        df_stream.emit(_df.compute())
        batch_progress_bar.update(1)
        #if partition > 2:
        #    break

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                sleep(10)
                break
        client.close()
        cluster.close()

    num_pairs_written = writer.close()
    batch_progress_bar.close()
    sys.stderr.write("\n\n")
    return num_pairs_written


def to_hic(df: AlignDf, max_fragment_id: int = 0) -> HicTxtDf:
    res = []

    keep_segments = (
        df.query("pass_filter == True").replace({"strand": {True: "0", False: "16"}})
        .astype({"strand": HicTxtDf.DTYPE['strand1'], "chrom": str})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):

        if max_fragment_id:
            mask = read_df.fragment_id > max_fragment_id
            if mask.any():
                logger.warning("Overflow found in fragment id, removing {} records" % mask.sum())
                read_df = read_df[~mask]

        if len(read_df) <= 1:
            continue
        read_id = "{}:{}".format(read_df.read_name.values[0], read_idx)
        rows = list(read_df.sort_values(["read_start"], ascending=True).itertuples())

        for (pos1, pos2) in combinations(rows, 2):
            if pos1.align_idx == pos2.align_idx:
                raise ValueError
            if pos1.fragment_id > pos2.fragment_id:
                pos1, pos2 = pos2, pos1
            res.append(
                (read_id,
                 pos1.strand,
                 pos1.chrom,
                 pos1.fragment_end,
                 pos1.fragment_id - 1,
                 pos2.strand,
                 pos2.chrom,
                 pos2.fragment_end,
                 pos2.fragment_id - 1,
                 pos1.mapping_quality,
                 pos2.mapping_quality,
                 )
            )
    return pd.DataFrame(res, columns=HicTxtDf.DTYPE.keys()).astype(HicTxtDf.DTYPE)

