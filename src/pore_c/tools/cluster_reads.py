import random
import sys
from typing import Optional, Tuple

import networkx as nx
import numpy as np
from intervaltree import IntervalTree
from pysam import AlignmentFile


def read_mappings_iter(bam, sort_flag=False, unique_intervals=False):
    try:
        assert sort_flag in ["start", "end", False]
    except:
        raise ValueError('sort flag must be either "start" or "end" or False.')
    aligns = []
    current_read_name = None
    seen = {}
    for align in bam:
        if align.is_unmapped:
            yield ([align])
            continue
        if current_read_name is None:
            current_read_name = align.query_name
            aligns.append(align)
        elif current_read_name == align.query_name:
            aligns.append(align)
        else:
            if unique_intervals:
                for entry in aligns:
                    st, en = entry.query_alignment_start, entry.query_alignment_end
                    if (st, en) not in seen:
                        seen[(st, en)] = entry
                    elif entry.get_tag("AS") > seen[(st, en)].get_tag("AS"):
                        seen[(st, en)] = entry
                aligns = seen.values()
                seen = {}
            if sort_flag == "end":
                yield sorted(aligns, key=lambda x: x.query_alignment_end)
            elif sort_flag == "start":
                yield sorted(aligns, key=lambda x: x.query_alignment_start)
            else:
                yield aligns
            current_read_name = align.query_name
            aligns = [align]
    if unique_intervals:
        for entry in aligns:
            st, en = entry.query_alignment_start, entry.query_alignment_end
            if (st, en) not in seen:
                seen[(st, en)] = entry
            elif entry.get_tag("AS") > seen[(st, en)].get_tag("AS"):
                seen[(st, en)] = entry
        aligns = seen.values()
    if sort_flag == "end":
        yield sorted(aligns, key=lambda x: x.query_alignment_end)
    elif sort_flag == "start":
        yield sorted(aligns, key=lambda x: x.query_alignment_start)
    else:
        yield aligns


def read_midpoint(align):
    return int((align.query_alignment_end + align.query_alignment_start) / 2)


def contained_segments_filter(aligns, mapping_quality_cutoff=0):
    """ takes a list of alignments (as gathered by the namesorted bam iterator).
    returns indices into that list of aligns to keep.
    """

    discard = set()

    intervals = [
        (align.query_alignment_start, align.query_alignment_end, x)
        for x, align in enumerate(aligns)
        if align.mapping_quality >= mapping_quality_cutoff
    ]

    tree = IntervalTree.from_tuples(intervals)

    for interval in tree:
        start, stop, val = interval
        for hits in tree.envelop(start + 1, stop - 1):
            new_discards = set(list(zip(hits))[2])
            discard = discard.union(new_discards)
    keep = set(range(len(aligns))) - discard
    return sorted(list(keep))


def cluster_aligned_segments(aligns, trim, mapping_quality_cutoff=0):
    if len(aligns) == 1:
        return []

    intervals = [
        (
            min(read_midpoint(align), align.query_alignment_start + trim),
            max(read_midpoint(align) + 1, align.query_alignment_end - trim),
            x,
        )
        for x, align in enumerate(aligns)
        if align.mapping_quality >= mapping_quality_cutoff
    ]

    tree = IntervalTree.from_tuples(intervals)
    tree.merge_overlaps(data_initializer=[], data_reducer=lambda x, y: x + [y])

    if len(tree) == 1:
        return []

    keep = []
    for interval in tree:
        indices = interval.data
        scores = np.array([aligns[x].get_tag("AS") for x in indices])
        # TODO: implement a threshold between best and second best scores to remove ambiguous mappings
        best_read = np.argmax(scores)
        keep.append(indices[best_read])
    return sorted(keep)


def overlap(align1, align2):
    start1 = align1.query_alignment_start
    end1 = align1.query_alignment_end
    start2 = align2.query_alignment_start
    end2 = align2.query_alignment_end

    if end1 < start2 or end2 < start1:
        return 0
    elif start1 <= start2 and end1 >= end2:
        return align1.query_alignment_length
    elif start1 >= start2 and end1 <= end2:
        return align2.query_alignment_length
    elif start1 <= start2 and end1 <= end2:
        return end1 - start2
    elif start1 >= start2 and end1 >= end2:
        return end2 - start1


def measure_overlaps(input_bam: str, output_table: str, no_zero):

    bam_in = AlignmentFile(input_bam)
    f_out = open(output_table, "w")
    header = "read_id,align1,align2,overlap_length\n"
    f_out.write(header)

    for read_aligns in read_mappings_iter(bam_in):
        l = len(read_aligns)
        for x in range(l - 1):
            for y in range(1, l):
                olap = overlap(read_aligns[x], read_aligns[y])
                if no_zero and olap == 0:
                    continue
                f_out.write(
                    "{read_id},{align1},{align2},{overlap}\n".format(
                        read_id=read_aligns[x].query_name, align1=x, align2=y, overlap=olap
                    )
                )


def remove_contained_segments(
    input_bam: str,
    keep_bam: str,
    discard_bam: str,
    mapping_quality_cutoff: int,
    alignment_stats: Optional[str] = None,
) -> Tuple[int, int, int, int]:

    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, "wb", template=bam_in)
    bam_discard = AlignmentFile(discard_bam, "wb", template=bam_in)

    if alignment_stats != None:
        alignment_stats_out = open(alignment_stats, "w")
        alignment_stats_out.write("read_id,mapping_id,filter_retained,query_start,query_end,mapq\n")

    for read_aligns in read_mappings_iter(bam_in):
        keep = contained_segments_filter(read_aligns, mapping_quality_cutoff=mapping_quality_cutoff)

        if len(keep) == 0:
            for idx, align in enumerate(read_aligns):
                bam_discard.write(align)
                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )
            continue

        for idx, align in enumerate(read_aligns):
            if idx in keep:
                bam_keep.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=1,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )
            else:
                bam_discard.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )


def cluster_reads(
    input_bam: str,
    keep_bam: str,
    discard_bam: str,
    trim: int,
    mapping_quality_cutoff: int,
    alignment_stats: Optional[str] = None,
) -> Tuple[int, int, int, int]:

    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, "wb", template=bam_in)
    bam_discard = AlignmentFile(discard_bam, "wb", template=bam_in)

    if alignment_stats != None:
        alignment_stats_out = open(alignment_stats, "w")
        alignment_stats_out.write("read_id,mapping_id,filter_retained,query_start,query_end,mapq\n")

    num_reads = 0
    num_reads_kept = 0
    num_aligns = 0
    num_aligns_kept = 0
    for read_aligns in read_mappings_iter(bam_in):
        # read_aligns is a list of pysam.AlignmentFile objects
        num_reads += 1
        num_aligns += len(read_aligns)

        # non_substring_aligns = [aligns[x] for x in substring_filter(read_aligns,mapping_quality_cutoff = mapping_quality_cutoff)]

        keep = cluster_aligned_segments(read_aligns, trim, mapping_quality_cutoff)
        # keep is a list of indices into the read_aligns object of which alignments to keep

        if len(keep) == 0:
            for idx, align in enumerate(read_aligns):
                bam_discard.write(align)
                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )

            continue

        for idx, align in enumerate(read_aligns):
            if idx in keep:
                bam_keep.write(read_aligns[idx])

                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=1,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )
            else:
                bam_discard.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )

        num_aligns_kept += len(keep)
    return (num_reads, num_reads_kept, num_aligns, num_aligns_kept)


##################
###fragdag code###
##################


def minimap_gapscore(length, o1=4, o2=24, e1=2, e2=1):
    return min([int(o1) + int(length) * int(e1), int(o2) + int(length) * int(e2)])


# O=5 E=2 is default for bwa bwasw
# O=6 E=1 is default for bwa mem


def bwa_gapscore(length, O=5, E=2):
    return O + length * E


def fragDAG(aligns, aligner="minimap2", params="default"):
    if len(aligns) == 0:
        return [], []
    G = nx.DiGraph()
    edge_values = {}
    G.add_node("IN")
    G.add_node("OUT")
    L = len(aligns)
    readlen = aligns[0].query_length
    for x in range(L):
        G.add_node(x)
        if aligner == "minimap2":
            GP_in = minimap_gapscore(aligns[x].query_alignment_start)
            GP_out = minimap_gapscore(readlen - aligns[x].query_alignment_end)
        elif aligner == "bwa":
            GP_in = bwa_gapscore(aligns[x].query_alignment_start)
            GP_out = bwa_gapscore(readlen - aligns[x].query_alignment_end)
        AS = aligns[x].get_tag("AS")
        edge_values[("IN", x)] = (GP_in, AS)  # store this data for diagnosis
        edge_values[(x, "OUT")] = (GP_out, AS)  # store this data for diagnosis
        G.add_edge("IN", x, weight=GP_in - AS)
        G.add_edge(x, "OUT", weight=GP_out)
    for x in range(L - 1):
        for y in range(x + 1, L):

            if aligner == "minimap2":
                GP = minimap_gapscore(
                    abs(aligns[y].query_alignment_start - aligns[x].query_alignment_end)
                )
            elif aligner == "bwa":
                GP = bwa_gapscore(
                    abs(aligns[y].query_alignment_start - aligns[x].query_alignment_end)
                )
            AS = aligns[y].get_tag("AS")
            edge_values[(x, y)] = (GP, AS)  # store this data for diagnosis
            G.add_edge(x, y, weight=GP - AS)

    # returns a list of the path through the graph as well as the
    return (
        nx.single_source_bellman_ford(G, "IN", "OUT")[1],
        edge_values,
    )  # including this as output enables reconstruction of the DAG that ended up filtering reads out. useful for diagnostic purposes


def fragDAG_filter(
    input_bam: str,
    keep_bam: str,
    discard_bam: str,
    mapping_quality_cutoff: int,
    aligner: str,
    aligner_params: Optional[str] = None,
    stats: Optional[str] = None,
    graph: Optional[str] = None,
):
    # the graph structure for each read should maybe be stored for diagnostics, but it isn't clear how to go about doing that tersely. maybe a string of the data dictionary?

    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, "wb", template=bam_in)
    bam_discard = AlignmentFile(discard_bam, "wb", template=bam_in)

    if stats != None:
        alignment_stats_out = open(stats, "w")
        alignment_stats_out.write("read_id,mapping_id,filter_retained,query_start,query_end,mapq\n")

    if graph != None:
        graph_stats_out = open(graph, "w")

    for _read_aligns in read_mappings_iter(bam_in, sort_flag="end", unique_intervals=True):

        # address unmapped reads
        if _read_aligns[0].is_unmapped:
            if stats != None:
                align = _read_aligns[0]
                alignment_stats_out.write(
                    "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                        read_id=align.query_name,
                        mapping_id=0,
                        filter_retained=-1,
                        q_start=0,
                        q_end=0,
                        mapq=align.mapq,
                    )
                )
            continue

        # filter on quality score, writing out anything that doesn't pass
        read_aligns = []
        for entry in _read_aligns:
            if entry.mapq < mapping_quality_cutoff:
                bam_discard.write(entry)
                if stats != None:
                    align = entry  # FIXME
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )
            else:
                read_aligns.append(entry)

        # this line returns two lists, the number of the reads that are part of the best scoring path, and the scores for that path based on the scoring function outlined
        if len(read_aligns) > 0:
            keep, graph_data = fragDAG(read_aligns, aligner=aligner, params=aligner_params)

        if len(keep) == 0:
            for idx, align in enumerate(read_aligns):
                bam_discard.write(align)
                if stats != None:
                    alignment_stats_out.write(
                        "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                            read_id=align.query_name,
                            mapping_id=idx,
                            filter_retained=0,
                            q_start=align.query_alignment_start,
                            q_end=align.query_alignment_end,
                            mapq=align.mapq,
                        )
                    )
            continue

        else:
            if graph != None:
                graph_stats_out.write("{}|{}\n".format(read_aligns[0].query_name, str(graph_data)))
            for idx, align in enumerate(read_aligns):
                if idx in keep:
                    bam_keep.write(read_aligns[idx])
                    if stats != None:
                        alignment_stats_out.write(
                            "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                                read_id=align.query_name,
                                mapping_id=idx,
                                filter_retained=1,
                                q_start=align.query_alignment_start,
                                q_end=align.query_alignment_end,
                                mapq=align.mapq,
                            )
                        )
                else:
                    bam_discard.write(read_aligns[idx])
                    if stats != None:
                        alignment_stats_out.write(
                            "{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n".format(
                                read_id=align.query_name,
                                mapping_id=idx,
                                filter_retained=0,
                                q_start=align.query_alignment_start,
                                q_end=align.query_alignment_end,
                                mapq=align.mapq,
                            )
                        )
