import logging
from itertools import combinations

import networkx as nx
import numpy as np

from pore_c.model import (
    AlignmentRecordDf,
    FragmentRecordDf,
    PoreCConcatemerRecord,
    PoreCConcatemerRecordDf,
    PoreCContactRecord,
    PoreCContactRecordDf,
    PoreCRecordDf,
)


logger = logging.getLogger(__name__)


def assign_fragments(
    align_table: AlignmentRecordDf,
    fragment_df: FragmentRecordDf,
    mapping_quality_cutoff: int = 1,
    min_overlap_length: int = 10,
    containment_cutoff: float = 99.0,
) -> PoreCRecordDf:
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
    # TODO: pyranges API has changed so that the new_position call overwrites the alignment start and end
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
    return align_df


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
        best_align_idx = read_df.loc[overlap_on_read, :].groupby(["read_start", "read_end"])["align_score"].idxmax()
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


def to_contacts(df: PoreCRecordDf) -> PoreCContactRecordDf:
    res = []
    keep_segments = df.query("pass_filter == True").assign(
        fragment_midpoint=lambda x: np.rint((x.fragment_start + x.fragment_end) * 0.5).astype(int)
    )
    chrom_dtype = df.dtypes["chrom"]
    for x, (read_idx, read_df) in enumerate(keep_segments.groupby("read_idx", as_index=False)):
        if len(read_df) <= 1:
            continue

        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .assign(pos_on_read=lambda x: np.arange(len(x)))
            .itertuples()
        )
        for (align_1, align_2) in combinations(rows, 2):
            contact = PoreCContactRecord.from_pore_c_align_pair(
                read_idx, align_1, align_2, contact_is_direct=align_2.pos_on_read - align_1.pos_on_read == 1
            )
            res.append(contact)

    res = PoreCContactRecord.to_dataframe(res, overrides={"align1_chrom": chrom_dtype, "align2_chrom": chrom_dtype})
    return res


def gather_concatemer_stats(contact_df: PoreCContactRecordDf) -> PoreCConcatemerRecordDf:

    by_read = contact_df.groupby(level="read_idx")
    num_reads = len(by_read)
    by_read_and_type = contact_df.reset_index().groupby(["read_idx", "contact_is_direct"])

    contact_counts = (
        by_read_and_type.size()
        .unstack(fill_value=0)
        .rename(columns={True: "direct_contacts", False: "indirect_contacts"})
    )
    contact_counts["total_contacts"] = contact_counts.sum(axis=1)
    contact_counts["read_order"] = contact_counts["direct_contacts"] + 1
    contact_counts["total_cis_contacts"] = by_read["contact_is_cis"].sum().astype(int)

    haplotype_stats = by_read["haplotype_pair_type"].value_counts().unstack(fill_value=0)
    drop = []
    for cat in contact_df.haplotype_pair_type.cat.categories:
        if cat == "null" or cat == "trans":
            if cat in haplotype_stats.columns:
                drop.append(cat)
        elif cat not in haplotype_stats.columns:
            haplotype_stats[cat] = 0
    if drop:
        haplotype_stats = haplotype_stats.drop(columns=drop)
    haplotype_stats = haplotype_stats.add_prefix("haplotype_")

    max_distance = (
        by_read_and_type[["contact_genome_distance", "contact_fragment_distance"]]
        .max()
        .unstack(fill_value=0)
        .astype(int)
        .rename(columns={True: "direct", False: "indirect"})
    )
    max_distance.columns = ["max_{1}_{0}".format(*_) for _ in max_distance.columns]
    dtype = PoreCConcatemerRecord.pandas_dtype()
    res = (
        contact_counts.join(haplotype_stats, how="left")
        .join(max_distance, how="left")
        .reset_index()
        .astype(dtype)[list(dtype.keys())]
    )
    assert len(res) == num_reads
    return res
