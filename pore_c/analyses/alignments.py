import logging
from itertools import combinations
from pathlib import Path
from typing import Dict, Union

import networkx as nx
import numpy as np
import pandas as pd

from pore_c.model import (
    AlignDf,
    BamEntryDf,
    ContactDf,
    FragmentDf,
    PairDf,
    PhasedBamEntryDf,
    PhasedContactDf,
    PhasedPoreCAlignDf,
    PoreCAlignDf,
    PoreCReadDf,
)
from pore_c.utils import DaskExecEnv, DataFrameProgress


logger = logging.getLogger(__name__)

FILTER_REASON_DTYPE = PoreCAlignDf.DTYPE["reason"]


def assign_fragments(
    align_table: Union[BamEntryDf, PhasedBamEntryDf],
    fragment_df: FragmentDf,
    mapping_quality_cutoff: int = 1,
    min_overlap_length: int = 10,
    containment_cutoff: float = 99.0,
) -> Dict:
    # initialise everything as having passed filter

    from pore_c.model import PoreCRecord

    pore_c_table = PoreCRecord.init_dataframe(align_table.reset_index(drop=False))
    dtype = pore_c_table.dtypes

    align_types = pore_c_table.align_type.value_counts()
    some_aligns = align_types["unmapped"] != align_types.sum()

    if not some_aligns:
        logger.warning("No alignments in dataframe")
        return pore_c_table

    fragment_assignments = assign_fragment(pore_c_table, fragment_df, min_overlap_length, containment_cutoff)

    pore_c_table = pore_c_table.set_index("align_idx", drop=False)
    pore_c_table.update(fragment_assignments)
    pore_c_table = pore_c_table.astype(dtype)

    # apply alignment-level filters
    unmapped_mask = pore_c_table.align_type == "unmapped"
    pore_c_table.loc[unmapped_mask, "pass_filter"] = False
    pore_c_table.loc[unmapped_mask, "filter_reason"] = "unmapped"

    # if not unmapped, but mapping quality below cutoff then fail
    fail_mq_mask = ~unmapped_mask & (pore_c_table.mapping_quality <= mapping_quality_cutoff)
    pore_c_table.loc[fail_mq_mask, "pass_filter"] = False
    pore_c_table.loc[fail_mq_mask, "filter_reason"] = "low_mq"

    # no need to do other checks if nothing left
    if pore_c_table["pass_filter"].any():
        # for the remaining alignments filtering happens on a per-read basis
        by_read_res = (
            pore_c_table[pore_c_table.pass_filter]
            .groupby("read_idx", sort=False, as_index=False)
            .apply(apply_per_read_filters)
        )
        pore_c_table.update(by_read_res[["pass_filter", "filter_reason"]])
        pore_c_table = pore_c_table.astype(dtype)
    else:
        logger.warning("No alignments passed filter")

    return pore_c_table.reset_index(drop=True).sort_values(["read_idx", "align_idx"], ascending=True)


def assign_fragment(pore_c_table, fragment_df, min_overlap_length: int, containment_cutoff: float):
    import pyranges as pr

    align_range = pr.PyRanges(
        pore_c_table[["chrom", "start", "end", "align_idx"]].rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )
    fragment_range = pr.PyRanges(
        fragment_df[["chrom", "start", "end", "fragment_id"]].rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )
    # all overlaps, one to many
    overlaps = align_range.join(fragment_range).new_position("intersection")
    if len(overlaps) == 0:
        raise ValueError("No overlaps found between alignments and fragments, this shouldn't happen")
    overlaps = (
        overlaps.df.rename(
            columns={
                "Start": "start",
                "End": "end",
                "Start_a": "align_start",
                "End_a": "align_end",
                "Start_b": "fragment_start",
                "End_b": "fragment_end",
            }
        )
        .eval("overlap_length = (end - start)")
        .query(f"overlap_length >= {min_overlap_length}")  # TODO: what if restriction fragment < minimum
        .eval("perc_of_alignment = (100.0 * overlap_length) / (align_end - align_start)")
        .eval("perc_of_fragment = (100.0 * overlap_length) / (fragment_end - fragment_start)")
        .eval(f"is_contained = (perc_of_fragment >= {containment_cutoff})")
    )

    # per-alignment statistics
    by_align = overlaps.groupby("align_idx", sort=True)

    rank = by_align["overlap_length"].rank(method="first", ascending=False).astype(int)
    overlaps["overlap_length_rank"] = rank

    best_overlap = overlaps[overlaps.overlap_length_rank == 1].set_index(["align_idx"])

    contained_fragments = (
        by_align["is_contained"]
        .agg(["size", "sum"])
        .astype({"sum": int})
        .rename(columns={"size": "num_overlapping_fragments", "sum": "num_contained_fragments"})
    )
    align_df = contained_fragments.join(
        best_overlap[
            [
                "fragment_id",
                "fragment_start",
                "fragment_end",
                "overlap_length",
                "perc_of_alignment",
                "perc_of_fragment",
                "is_contained",
            ]
        ]
    )

    dtype = {col: dtype for col, dtype in pore_c_table.dtypes.items() if col in align_df.columns}
    align_df = align_df.astype(dtype)
    # pore_c_table.head(5))
    return align_df


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
    return read_df.pipe(filter_singleton).pipe(filter_exact_overlap_on_query).pipe(filter_shortest_path)


def filter_singleton(read_df):
    if len(read_df) == 1:  # if you have a single alignment at this point you fail
        read_df.loc[:, "pass_filter"] = False
        read_df.loc[:, "filter_reason"] = "singleton"
    return read_df


def filter_exact_overlap_on_query(read_df):
    overlap_on_read = read_df.duplicated(subset=["read_start", "read_end"], keep=False)
    if overlap_on_read.any():
        best_align_idx = read_df.loc[overlap_on_read, :].groupby(["read_start", "read_end"])["score"].idxmax()
        overlap_on_read[best_align_idx.values] = False
        read_df.loc[overlap_on_read, "pass_filter"] = False
        read_df.loc[overlap_on_read, "filter_reason"] = "overlap_on_read"
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
    aligns = aligns[["read_start", "read_end", "read_length", "align_score"]].copy().sort_values(["read_end"])
    node_ids = list(aligns.index)
    graph = nx.DiGraph()
    # initialise graph with root and sink node, and one for each alignment
    # edges in the graph will represent transitions from one alignment segment
    # to the next
    graph.add_nodes_from(["ROOT", "SINK"] + node_ids)
    for align in aligns.itertuples():
        align_idx = align.Index
        align_score = align.align_score
        gap_penalty_start = gap_fn(align.read_start)
        graph.add_edge("ROOT", align_idx, weight=gap_penalty_start - align_score)
        gap_penalty_end = gap_fn(int(align.read_length - align.read_end))
        graph.add_edge(align_idx, "SINK", weight=gap_penalty_end)

    # for each pair of aligned segments add an edge
    for idx_a, align_idx_a in enumerate(node_ids[:-1]):
        align_a_end = aligns.at[align_idx_a, "read_end"]
        for align_idx_b in node_ids[idx_a + 1 :]:  # noqa: E203 black does this
            align_b_score = aligns.at[align_idx_b, "align_score"]
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
            read_df.at[idx, "filter_reason"] = "not_on_shortest_path"
    return read_df


def clean_haplotypes(align_df, output_path, n_workers: int = 1):

    with DaskExecEnv(n_workers=n_workers):
        cleaned_df = align_df.map_partitions(_clean_haplotypes, meta=align_df._meta).compute()
    return cleaned_df


def _clean_haplotypes(align_df, *args):

    return align_df.groupby("read_idx").apply(clean_read_haplotypes)


def clean_read_haplotypes(read_aligns, min_count=1, min_prop=0.5):

    # dataframe containing the count of phased alignment segments per chromosome and haplotype
    counts_df = (
        read_aligns.loc[(read_aligns.pass_filter & (read_aligns.haplotype > 0)), :]
        .groupby(["chrom", "phase_block", "haplotype"])
        .size()
        .rename("count")
        .to_frame()
        .groupby(level=["chrom", "phase_block"])
        .filter(lambda x: len(x) > 1)
    )
    if len(counts_df) == 0:
        return read_aligns
    else:
        g = counts_df.groupby(level=["chrom", "phase_block"])
        # proportion of alignments in each phase block assigned to each haplotype
        counts_df["prop"] = g.transform(lambda x: x / x.sum())
        counts_df["rank"] = g["count"].rank()
        raise NotImplementedError
        return read_aligns


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


def convert_align_df_to_contact_df(
    align_df: Union[PoreCAlignDf, PhasedPoreCAlignDf], contacts: Path, n_workers: int = 1
):

    phased = "phase_block" in align_df.columns
    if phased:
        df_type = PhasedContactDf
    else:
        df_type = ContactDf
    with DaskExecEnv(n_workers=n_workers):
        contact_df = align_df.map_partitions(to_contacts, phased=phased, meta=df_type.DTYPE)
        contact_df.to_parquet(str(contacts), engine="pyarrow")


def to_contacts(df: Union[AlignDf, PhasedPoreCAlignDf], phased=False) -> ContactDf:

    res = []
    keep_segments = (
        df.query("pass_filter == True").replace({"strand": {True: "+", False: "-"}})
        # convert strand to +/- and chrom to string for lexographical sorting to get upper-triangle
        .astype({"strand": PairDf.DTYPE["strand1"], "chrom": str})
    )
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue

        read_info = read_df[["read_name", "read_length"]].iloc[0, :]
        read_name = read_info["read_name"]
        read_length = read_info["read_length"]
        read_order = len(read_df)
        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .assign(
                pos_on_read=lambda x: np.arange(len(x)),
                fragment_midpoint=lambda x: np.rint((x.fragment_start + x.fragment_end) * 0.5).astype(int),
            )
            .itertuples()
        )

        for (align_1, align_2) in combinations(rows, 2):
            res.append(
                align_pair_to_tuple(align_1, align_2, read_name, read_idx, read_length, read_order, phased=phased)
            )
    if phased:
        df_type = PhasedContactDf
    else:
        df_type = ContactDf
    return pd.DataFrame(res, columns=df_type.DTYPE.keys()).astype(df_type.DTYPE)


def align_pair_to_tuple(align_1, align_2, read_name, read_idx, read_length, read_order, phased=False):
    contact_is_direct = (align_2.pos_on_read - align_1.pos_on_read) == 1
    if align_1.fragment_id > align_2.fragment_id:
        align_1, align_2 = align_2, align_1

    contact_is_cis = align_1.chrom == align_2.chrom
    if contact_is_cis:
        contact_genome_distance = align_2.start - align_1.end
        contact_fragment_distance = align_2.fragment_midpoint - align_1.fragment_midpoint
    else:
        contact_genome_distance = 0
        contact_fragment_distance = 0

    res = (
        read_name,
        read_idx,
        read_length,
        read_order,
        contact_is_direct,
        contact_is_cis,
        contact_genome_distance,
        contact_fragment_distance,
        align_1.align_idx,
        align_1.chrom,
        align_1.start,
        align_1.end,
        align_1.strand,
        align_1.read_start,
        align_1.read_end,
        align_1.mapping_quality,
        align_1.score,
        align_1.fragment_id,
        align_1.fragment_start,
        align_1.fragment_end,
        align_1.fragment_midpoint,
        align_2.align_idx,
        align_2.chrom,
        align_2.start,
        align_2.end,
        align_2.strand,
        align_2.read_start,
        align_2.read_end,
        align_2.mapping_quality,
        align_2.score,
        align_2.fragment_id,
        align_2.fragment_start,
        align_2.fragment_end,
        align_2.fragment_midpoint,
    )
    if phased:
        res = (*res, align_1.phase_block, align_1.haplotype, align_2.phase_block, align_2.haplotype)
    return res
