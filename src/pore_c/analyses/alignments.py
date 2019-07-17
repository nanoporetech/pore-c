import logging
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from streamz import Stream
from tqdm import tqdm

from pore_c.datasources import NameSortedBamSource
from pore_c.io import TableWriter
from pore_c.model import (BamEntryDf, FragmentDf,
                          PoreCAlignDf, PoreCReadDf)
from pore_c.utils import DataFrameProgress

logger = logging.getLogger(__name__)

FILTER_REASON_DTYPE = PoreCAlignDf.DTYPE["reason"]


def parse_alignment_bam(
    input_bam: Path,
    fragment_df: FragmentDf,
    alignment_table: Path = None,
    read_table: Path = None,
    overlap_table: Path = None,
    alignment_summary: Path = None,
    read_summary: Path = None,
    chunksize: int = 50000,
    n_workers: int = 1,
):
    """Filter alignments to keep only alignments that contribute to contacts

    Parameters
    ----------

    input_bam : str
                Path to a namesorted bam with unfiltered alignments
    chunksize: int
                The alignments are batched for processing, this controls the batch size

    """

    source_aligns = NameSortedBamSource(input_bam, metadata={})
    source_aligns.discover()

    parallel = n_workers > 1
    fragment_df = fragment_df.set_index(["fragment_id"]).sort_index()  # .rename_axis("index", axis=0)
    if parallel:
        from dask.distributed import Client, LocalCluster
        from time import sleep

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)
        fragment_df = client.scatter(fragment_df)

    writers = dict(
        alignment_table=TableWriter(alignment_table),
        read_table=TableWriter(read_table),
        overlap_table=TableWriter(overlap_table),
    )

    batch_progress_bar = tqdm(total=None, desc="Alignments submitted: ", unit=" alignments", position=0)
    alignment_progress = AlignmentProgress(position=1)
    read_progress = ReadProgress(position=2)
    # perc_alignment_bar = tqdm(total=None, desc="Alignments processed: ", unit=" alignments", position=1)

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
    align_sink = (  # noqa: F841
        filtered_align_stream.pluck("alignment_table")
        .accumulate(alignment_progress, returns_state=True, start=alignment_progress)
        .sink(writers["alignment_table"])
    )

    read_sink = (  # noqa: F841
        filtered_align_stream.pluck("read_table")
        .accumulate(read_progress, returns_state=True, start=read_progress)
        .sink(writers["read_table"])
    )

    overlap_sink = filtered_align_stream.pluck("overlap_table").sink(writers["overlap_table"])  # noqa: F841

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
    alignment_progress.close()
    alignment_progress.save(alignment_summary)
    read_progress.close()
    read_progress.save(read_summary)
    sys.stderr.write("\n\n\n")
    sys.stdout.write("\n")
    return read_progress.final_stats()


def filter_read_alignments(
    df: BamEntryDf,
    fragment_df: FragmentDf,
    mapping_quality_cutoff: int = 1,
    min_overlap_length: int = 10,
    containment_cutoff: float = 99.0,
):
    # initialise everything as having passed filter
    res = df.assign(pass_filter=True, reason="pass").astype({"reason": FILTER_REASON_DTYPE}).set_index("align_idx")

    # calculate fragment overlaps
    overlaps = (
        res.ginterval.overlap(fragment_df)
        .query("overlap_length >= {}".format(min_overlap_length))
        .loc[
            :,
            ["align_idx", "fragment_id", "other_start", "other_end", "overlap_length", "perc_of_self", "perc_of_other"],
        ]
        .rename(
            columns={
                "perc_of_self": "perc_of_alignment",
                "perc_of_other": "perc_of_fragment",
                "other_start": "fragment_start",
                "other_end": "fragment_end",
            }
        )
        .sort_values(["overlap_length"], ascending=False)
    )
    overlaps = pd.merge(overlaps, res[["read_idx"]], left_on="align_idx", right_index=True, how="left")

    contained_fragments = (
        overlaps.query("perc_of_fragment >= {}".format(containment_cutoff))
        .groupby("align_idx", as_index=True, sort=False)
        .size()
        .rename("contained_fragments")
    )

    # pick the longest overlap as the associated fragment
    olg = overlaps.groupby("align_idx", as_index=True, sort=False)
    overlap_summary = (
        olg.head(1).set_index(["align_idx"]).join(contained_fragments)  # take longest overlap as fragment assignment
    )
    res = res.join(
        overlap_summary[
            [
                "fragment_id",
                "contained_fragments",
                "fragment_start",
                "fragment_end",
                "perc_of_alignment",
                "perc_of_fragment",
            ]
        ]
    ).fillna({"fragment_id": -1, "contained_fragments": 0, "fragment_start": 0, "fragment_end": 0})

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
            res.query("pass_filter == True")
            .groupby("read_name", sort=False, as_index=False)
            .apply(apply_per_read_filters)
        )
        res.update(by_read_res[["pass_filter", "reason"]])
        res = res.astype({"reason": FILTER_REASON_DTYPE})
    res = res.reset_index().porec_align.cast(subset=True, fillna=True)
    read_table = calculate_read_stats(res)

    return {
        "alignment_table": res.reset_index(),
        "overlap_table": overlaps.reset_index(),
        "read_table": read_table.reset_index(),
    }


def calculate_read_stats(align_df: PoreCAlignDf) -> PoreCReadDf:
    # total number of reads and how many aligments per read
    num_aligns = (
        align_df.groupby(["read_idx"], sort=False)[["read_name", "read_length"]]
        .max()
        .join(align_df.query("chrom != 'NULL'").groupby(["read_idx"]).size().rename("num_aligns").to_frame())
        .fillna(0)
        .astype({"num_aligns": np.uint8})
    )

    # subset of alignments that passed filters
    pass_aligns = align_df.query("pass_filter == True").eval(
        "perc_read_assigned = 100.0 * (read_end - read_start) / read_length"
    )

    pass_stats = (
        pass_aligns.groupby(["read_idx"], sort=False)
        .agg(
            {"fragment_id": "nunique", "pass_filter": "size", "contained_fragments": "sum", "perc_read_assigned": "sum"}
        )
        .rename(columns={"fragment_id": "unique_fragments_assigned", "pass_filter": "num_pass_aligns"})
    )

    def count_contacts(chroms: pd.Series):
        """Take a list of chromosomes and calculate the total number of contacts (N choose 2) and the
        number of those that are cis
        """
        if len(chroms) < 2:
            res = {"num_contacts": 0, "num_cis_contacts": 0}
        else:
            is_cis = np.array([chrom1 == chrom2 for (chrom1, chrom2) in combinations(chroms, 2)])
            res = {"num_contacts": len(is_cis), "num_cis_contacts": is_cis.sum()}
        res["num_chroms_contacted"] = chroms.nunique()
        return pd.DataFrame(res, index=[chroms.name])

    def _check_index(_df):
        level_0 = _df.index.get_level_values(0)
        level_1 = _df.index.get_level_values(1)
        assert (level_0 == level_1).all()
        return _df

    contact_stats = (
        pass_aligns.groupby(["read_idx"], sort=False)["chrom"].apply(count_contacts).pipe(_check_index).droplevel(1)
    )

    # create a merged dataframe with one row per read
    res = num_aligns.join(pass_stats.join(contact_stats), how="outer")
    return res.reset_index().porec_read.cast(subset=True, fillna=True)


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


def bwa_gapscore(length, O=5, E=2):  # noqa: E741
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
        for align_idx_b in node_ids[idx_a + 1:]:
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


class AlignmentProgress(DataFrameProgress):
    def __init__(self, **kwds):
        kwds["desc"] = "Alignments processed"
        kwds["unit"] = " alignments"
        super(AlignmentProgress, self).__init__(**kwds)

    def update_data(self, align_df):
        pass_stats = align_df["reason"].value_counts(dropna=False).to_frame().T.reset_index(drop=True)
        if self._data is None:
            self._data = pass_stats
        else:
            self._data += pass_stats

    def update_progress_bar(self, num_aligns):
        percents = (100.0 * self._data).div(self._data.sum(axis=1), axis=0)
        self._bar.update(num_aligns)
        self._bar.set_postfix(percents.iloc[0, :].to_dict())


class ReadProgress(DataFrameProgress):
    def __init__(self, **kwds):
        kwds["desc"] = "Reads processed"
        kwds["unit"] = " reads"
        super(ReadProgress, self).__init__(**kwds)
        self._contact_histogram = None

    def update_data(self, read_df):
        count_cols = ["read_length", "num_contacts", "num_cis_contacts", "num_aligns", "num_pass_aligns"]
        counts = read_df.loc[:, count_cols].sum(axis=0)
        counts["reads"] = len(read_df)
        counts = counts.astype(int).to_frame().T

        if self._data is None:
            self._data = counts
        else:
            self._data += counts

    @staticmethod
    def summarize(df):
        summary = (
            df.eval("Gb = read_length * 1e-9")
            .eval("contacts_per_Gb = num_contacts/Gb")
            .eval("percent_cis = 100 * (num_cis_contacts / num_contacts)")
            .loc[:, ["reads", "Gb", "num_contacts", "contacts_per_Gb", "percent_cis"]]
        )
        return summary

    def update_progress_bar(self, num_aligns):
        self._bar.update(num_aligns)
        self._bar.set_postfix(self.summarize(self._data).iloc[0, :].to_dict())

    def final_stats(self):
        return self._data.join(self.summarize(self._data).drop(columns=["num_contacts", "reads"])).iloc[0, :].to_dict()

    def save(self, path):
        df = self._data.join(self.summarize(self._data).drop(columns=["num_contacts", "reads"]))
        df.to_csv(path, index=False)
