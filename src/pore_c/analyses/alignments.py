import logging
from pathlib import Path
from typing import Dict, NewType, Union

import intake
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow as pa
import yaml
from intake.catalog.local import YAMLFileCatalog
from pyarrow import parquet as pq
from pysam import AlignmentFile
from tqdm import tqdm
from itertools import combinations

from pore_c.datasources import NameSortedBamSource
from pore_c.model import GenomeIntervalDf

logger = logging.getLogger(__name__)

FILTER_REASON_DTYPE = pd.CategoricalDtype(
    ["pass", "unmapped", "singleton", "low_mq", "overlap_on_read", "not_on_shortest_path"], ordered=True
)


class TableWriter(object):
    def __init__(self, path):
        self.path = path
        self._writer = None
        self._schema = None

    def write(self, df):
        table = pa.Table.from_pandas(df, preserve_index=False, schema=self._schema)
        if self._writer is None:
            self._writer = pq.ParquetWriter(self.path, schema=table.schema)
            self._schema = table.schema
        try:
            self._writer.write_table(table)
        except:
            print(df)
            raise

    def __call__(self, *args, **kwds):
        return self.write(*args, **kwds)

    def close(self):
        self._writer.close()


def filter_alignments(
    input_bam: str,
    fragment_df: pd.DataFrame,
    align_table: str = None,
    read_table: str = None,
    overlap_table: str = None,
    catalog_file: str = None,
    chunksize=50000,
    n_workers=10,
):
    """Filter alignments to keep only alignments that contribute to contacts

    Parameters
    ----------

    input_bam : str
                Path to a namesorted bam with unfiltered alignments
    chunksize: int
                The alignments are batched for processing, this controls the batch size

    """
    from streamz import Stream
    from time import sleep
    from dask.distributed import Client, LocalCluster
    from scipy.special import binom

    parallel = n_workers > 1
    fragment_df = fragment_df.set_index(["fragment_id"]).sort_index()  # .rename_axis("index", axis=0)
    if parallel:
        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)
        fragment_df = client.scatter(fragment_df)

    writers = dict(
        align_table=TableWriter(align_table),
        read_table=TableWriter(read_table),
        overlap_table=TableWriter(overlap_table),
    )

    batch_progress_bar = tqdm(total=None, desc="Alignments submitted: ", unit=" alignments", position=0)
    perc_alignment_bar = tqdm(total=None, desc="Alignments processed: ", unit=" alignments", position=1)
    perc_read_bar = tqdm(total=None, desc="Reads processed: ", unit=" reads", position=2)

    def alignment_progress(state, align_df, progress_bar=None):
        pass_stats = align_df["reason"].value_counts(dropna=False)
        if state is None:
            state = pass_stats
        else:
            state += pass_stats
        counts = state.sort_index().to_frame().T
        percents = (100.0 * counts).div(counts.sum(axis=1), axis=0)
        if progress_bar is not None:
            progress_bar.update(len(align_df))
            progress_bar.set_postfix(percents.loc["reason", :].to_dict())
        return state, align_df

    def read_progress(state, read_df, progress_bar=None):
        have_contacts = read_df["num_pass_aligns"] >= 2
        batch_stats = pd.Series(
            {
                "reads": len(read_df),
                "Gb": read_df.read_length.sum() * 1e-9,
                "reads_with_contacts": have_contacts.sum(),
                #'pass_aligns': read_df['num_pass_aligns'].sum(),
                "contacts": binom(read_df.loc[have_contacts, "num_pass_aligns"].values, 2).sum(),
                #'contacts': binom(read_df.loc[have_contacts, 'num_pass_aligns'].values, 2).sum(),
            }
        )
        batch_stats["contacts_per_Gb"] = batch_stats["contacts"] / batch_stats["Gb"]
        if state is None:
            state = batch_stats
        else:
            # print(state, batch_stats)
            state += batch_stats
        if progress_bar is not None:
            progress_bar.update(len(read_df))
            progress_bar.set_postfix(state.to_dict())
        return state, read_df

    # stream that holds the raw alignment dfs
    bam_stream = Stream()

    # stream that holds the filtered/processed alignments
    if parallel:
        filtered_align_stream = (
            bam_stream.scatter().map(filter_read_alignments, fragment_df=fragment_df).buffer(n_workers).gather()
        )
    else:
        filtered_align_stream = bam_stream.map(filter_read_alignments, fragment_df=fragment_df)

    # write the alignments using the table writer, updating progress bar as we go
    align_sink = (
        filtered_align_stream.pluck("alignment_table")
        .accumulate(alignment_progress, returns_state=True, start=None, progress_bar=perc_alignment_bar)
        .sink(writers["align_table"])
    )

    read_sink = (
        filtered_align_stream.pluck("read_table")
        .accumulate(read_progress, returns_state=True, start=None, progress_bar=perc_read_bar)
        .sink(writers["read_table"])
    )

    overlap_sink = filtered_align_stream.pluck("overlap_table").sink(writers["overlap_table"])

    source_aligns = NameSortedBamSource(input_bam, metadata={})
    for batch_idx, align_df in enumerate(source_aligns.read_chunked(chunksize=chunksize)):
        bam_stream.emit(align_df)
        batch_progress_bar.update(len(align_df))
        batch_progress_bar.set_postfix({"batches": batch_idx})

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
    perc_alignment_bar.close()
    perc_read_bar.close()
    create_catalog_yaml(writers, catalog_file)

    print("\nCompleted\n")


def create_catalog_yaml(writers, catalog_file):
    catalog = {
        "name": "pore_c_filtered_files",
        "description": "Output files of pore-c tools alignment filtering",
        "sources": {},
    }
    for file_key, writer in writers.items():
        if not writer:
            continue
        entry = {"args": {}}
        writer.close()
        if isinstance(writer, TableWriter):
            entry["args"]["urlpath"] = "{{ CATALOG_DIR }}/" + str(writer.path.name)
            entry["driver"] = "parquet"
        elif isinstance(writer, AlignmentFile):
            entry["args"]["urlpath"] = "{{ CATALOG_DIR }}/" + Path(writer.filename.decode("utf8")).name
            entry["driver"] = "pore_c.datasources.NameSortedBamSource"
        catalog["sources"][file_key] = entry

    with open(catalog_file, "w") as fh:
        fh.write(yaml.dump(catalog))
    return YAMLFileCatalog(str(catalog_file))


def calculate_read_stats(read_df):

    # total number of reads and how many aligments per read
    num_aligns = (
        read_df.groupby(["read_name"])[["read_length"]]
        .max()
        .join(read_df.query("chrom != 'NULL'").groupby(["read_name"]).size().rename("num_aligns").to_frame())
        .fillna(0)
        .astype({"num_aligns": np.uint8})
    )

    # subset of alignments that passed filters
    pass_aligns = read_df.query("pass_filter == True").eval(
        "perc_read_aligned = 100.0 * (read_end - read_start) / read_length"
    )

    # how many alignments pass filtering and how many unique chromosomes
    # do they hit
    def prop_cis(chroms):
        if len(chroms) < 2:
            return np.NaN
        is_cis = np.array([chrom1 == chrom2 for (chrom1, chrom2) in combinations(chroms, 2)])
        return is_cis.mean()

    read_stats = (
        pass_aligns.groupby(["read_name"])["chrom"]
        .agg(["nunique", "size", prop_cis])
        .rename(columns={"nunique": "num_chroms", "size": "num_pass_aligns"})
    )

    # figure out what the "main" chromosome is and how much of the read hits there
    main_chrom = (
        pass_aligns.groupby(["read_name", "chrom"])
        .agg({"perc_read_aligned": "sum", "read_end": "size"})
        .sort_values("perc_read_aligned", ascending=False)
        .groupby("read_name")
        .head(1)
        .reset_index()
        .rename(
            columns={
                "chrom": "main_chrom",
                "perc_read_aligned": "main_chrom_perc_read_aligned",
                "read_end": "main_chrom_num_aligns",
            }
        )
        .set_index(["read_name"])
    )

    # create a merged dataframe with one row per read
    res = num_aligns.join(read_stats.join(main_chrom), how="outer")

    unmapped = res["main_chrom"].isnull()
    res.loc[unmapped, "main_chrom"] = "NULL"
    res.loc[unmapped, "num_chroms"] = 0
    res.loc[unmapped, "num_pass_aligns"] = 0
    res.loc[unmapped, "main_chrom_num_aligns"] = 0
    res.loc[unmapped, "main_chrom_perc_read_aligned"] = 0.0
    res = res.astype(
        {"num_aligns": np.uint8, "num_chroms": np.uint8, "num_pass_aligns": np.uint8, "main_chrom_num_aligns": np.uint8}
    )
    return res


def filter_read_alignments(df, mapping_quality_cutoff: int = 1, fragment_df: GenomeIntervalDf = None, min_overlap_length: int = 10, containment_cutoff: float = 99.0):
    # initialise everything as having passed filter
    res = df.assign(pass_filter=True, reason="pass").astype({"reason": FILTER_REASON_DTYPE}).set_index("align_idx")

    # calculate fragment overlaps
    overlaps = (
        res.ginterval.overlap(fragment_df)
        .query("overlap_length >= {}".format(min_overlap_length))
        .loc[
            :,
            [
                "align_idx",
                "fragment_id",
                "other_start",
                "other_end",
                "overlap_length",
                "perc_of_self",
                "perc_of_other",
            ],
        ]
        .rename(columns={"perc_of_self": "perc_of_alignment", "perc_of_other": "perc_of_fragment", "other_start": "fragment_start", "other_end": "fragment_end"})
        .sort_values(["overlap_length"], ascending=False)
    )
    overlaps = pd.merge(overlaps, res[['read_idx']], left_on="align_idx", right_index=True, how="left")

    contained_fragments = overlaps.query("perc_of_fragment >= {}".format(containment_cutoff)).groupby("align_idx", as_index=True, sort=False).size().rename("contained_fragments")

    # pick the longest overlap as the associated fragment
    olg = overlaps.groupby("align_idx", as_index=True, sort=False)
    overlap_summary = (
        olg.head(1)  # take longest overlap as fragment assignment
        .set_index(["align_idx"])
        .join(contained_fragments)
    )
    res = res.join(overlap_summary[['fragment_id', 'contained_fragments', 'fragment_start', 'fragment_end', "perc_of_alignment", "perc_of_fragment"]]).fillna({'fragment_id': -1, 'contained_fragments': 0, 'fragment_start': 0, 'fragment_end': 0})

    # apply alignment-level filters
    unmapped_mask = res.mapping_type == "unmapped"
    res.loc[unmapped_mask, "pass_filter"] = False
    res.loc[unmapped_mask, "reason"] = "unmapped"

    # if not unmapped, but mapping quality below cutoff then fail
    fail_mq_mask = ~unmapped_mask & (df.mapping_quality <= mapping_quality_cutoff)
    res.loc[fail_mq_mask, "pass_filter"] = False
    res.loc[fail_mq_mask, "reason"] = "low_mq"

    # no need to do other checks if nothing left
    if res["pass_filter"].any():
        # for the remaining alignments filtering happens on a per-read basis
        by_read_res = (
            res.query("pass_filter == True").groupby("read_name", sort=False, as_index=False).apply(apply_per_read_filters)
        )
        res.update(by_read_res[["pass_filter", "reason"]])
        res = res.astype({"reason": FILTER_REASON_DTYPE})
    read_table = calculate_read_stats(res)

    return {
        "alignment_table": res.reset_index(),
        "overlap_table": overlaps.reset_index(),
        "read_table": read_table.reset_index(),
    }


def apply_per_read_filters(read_df):
    return read_df.pipe(filter_singleton).pipe(filter_overlap_on_query).pipe(filter_shortest_path)


def filter_singleton(read_df):
    if len(read_df) == 1:  # if you have a single alignment at this point you fail
        read_df.loc[:, "pass_filter"] = False
        read_df.loc[:, "reason"] = "singleton"
    return read_df


def filter_overlap_on_query(read_df):
    overlap_on_read = read_df.duplicated(subset=["read_start", "read_end"], keep=False)
    if overlap_on_read.any():
        best_align_idx = read_df.loc[overlap_on_read, :].groupby(["read_start", "read_end"])["score"].idxmax()
        overlap_on_read[best_align_idx.values] = False
        read_df.loc[overlap_on_read, "pass_filter"] = False
        read_df.loc[overlap_on_read, "reason"] = "overlap_on_read"
    return read_df


def minimap_gapscore(length, o1=4, o2=24, e1=2, e2=1):
    return min([o1 + int(length) * e1, o2 + int(length) * e2])


def bwa_gapscore(length, O=5, E=2):
    # O=5 E=2 is default for bwa bwasw
    # O=6 E=1 is default for bwa mem
    return O + length * E


def create_align_graph(aligns, gap_fn):
    # we'll visit the alignments in order of increasing endpoint on the read, need to keep
    # the ids as the index in the original list of aligns for filtering later
    aligns = aligns[["read_start", "read_end", "read_length", "score"]].copy().sort_values(["read_end"])
    node_ids = list(aligns.index)
    graph = nx.DiGraph()
    # initialise graph with root and sink node, and one for each alignment
    # edges in the graph will represent transitions from one alignment segment
    # to the next
    graph.add_nodes_from(["ROOT", "SINK"] + node_ids)
    for align in aligns.itertuples():
        align_idx = align.Index
        align_score = align.score
        gap_penalty_start = gap_fn(align.read_start)
        graph.add_edge("ROOT", align_idx, weight=gap_penalty_start - align_score)
        gap_penalty_end = gap_fn(int(align.read_length - align.read_end))
        graph.add_edge(align_idx, "SINK", weight=gap_penalty_end)

    # for each pair of aligned segments add an edge
    for idx_a, align_idx_a in enumerate(node_ids[:-1]):
        align_a_end = aligns.at[align_idx_a, "read_end"]
        for align_idx_b in node_ids[idx_a + 1 :]:
            align_b_score = aligns.at[align_idx_b, "score"]
            align_b_read_start = aligns.at[align_idx_b, "read_start"]
            gap_penalty = gap_fn(abs(int(align_b_read_start) - int(align_a_end)))
            graph.add_edge(align_idx_a, align_idx_b, weight=gap_penalty - align_b_score)

    return graph


def filter_shortest_path(read_df, aligner="minimap2"):
    aligns = read_df[read_df.pass_filter]
    num_aligns = len(aligns)
    if num_aligns < 2:
        # can't build a graph, so nothing is filtered by this method
        return read_df

    if aligner == "minimap2":
        gap_fn = minimap_gapscore
    elif aligner == "bwa":
        gap_fn = bwa_gapscore
    else:
        raise ValueError(f"Unrecognised aligner: {aligner}")
    graph = create_align_graph(aligns, gap_fn)
    distance, shortest_path = nx.single_source_bellman_ford(graph, "ROOT", "SINK")
    for idx in aligns.index:
        if idx not in shortest_path:
            read_df.at[idx, "pass_filter"] = False
            read_df.at[idx, "reason"] = "not_on_shortest_path"
    return read_df
