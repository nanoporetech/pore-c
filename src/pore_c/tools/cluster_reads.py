import sys
from typing import Tuple, Optional
import numpy as np
from pysam import AlignmentFile
from intervaltree import IntervalTree

import random

def read_mappings_iter(bam, mapping_quality_cutoff=0):
    aligns = []
    current_read_name = None
    for align in bam:
        if align.is_unmapped:
            continue
        if current_read_name is None:
            current_read_name = align.query_name
            aligns.append(align)
        elif current_read_name == align.query_name:
            aligns.append(align)
        else:
            yield aligns
            current_read_name = align.query_name
            aligns = [align]
    yield aligns

def read_midpoint(align):
    return int ((align.query_alignment_end + align.query_alignment_start ) / 2)

def contained_segments_filter(aligns, mapping_quality_cutoff = 0):
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
        start,stop,val = interval
        for hits in tree.envelop(start+1,stop-1):
            new_discards = set(list(zip(hits))[2])
            discard = discard.union(new_discards)
    keep = set(range(len(aligns))) - discard
    return sorted(list(keep))


def cluster_aligned_segments(aligns, trim, mapping_quality_cutoff=0):
    if len(aligns) == 1:
        return []

    intervals = [

        (min( read_midpoint( align ) , align.query_alignment_start + trim), max(read_midpoint(align) + 1, align.query_alignment_end - trim), x)
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
        #TODO: implement a threshold between best and second best scores to remove ambiguous mappings
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
        return  align1.query_alignment_length
    elif start1 >= start2 and end1 <= end2:
        return  align2.query_alignment_length
    elif start1 <= start2 and end1 <= end2:
        return end1 - start2
    elif start1 >= start2 and end1 >= end2:
        return end2 - start1

def measure_overlaps(input_bam: str, output_table: str):

    bam_in = AlignmentFile(input_bam)
    f_out = open(output_table,'w')
    header = "read_id,align1,align2,overlap_length\n"
    f_out.write(header)

    for read_aligns in read_mappings_iter(bam_in):
        l = len(read_aligns)
        for x in range(l-1):
            for y in range(1,l):
                f_out.write('{read_id},{align1},{align2},{overlap}\n'.format(read_id = read_aligns[x].query_name,align1 = x,align2 = y, overlap = overlap(read_aligns[x],read_aligns[y])))


def remove_contained_segments(input_bam: str, keep_bam: str, discard_bam: str, mapping_quality_cutoff: int, alignment_stats: Optional[str] = None) -> Tuple[int, int, int, int]:

    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, 'wb', template=bam_in)
    bam_discard = AlignmentFile(discard_bam, 'wb', template=bam_in)

    if alignment_stats != None:
        alignment_stats_out = open(alignment_stats,'w')
        alignment_stats_out.write('read_id,mapping_id,filter_retained,query_start,query_end,mapq\n')


    for read_aligns in read_mappings_iter(bam_in):
        keep = contained_segments_filter(read_aligns,mapping_quality_cutoff = mapping_quality_cutoff)

        if len(keep) == 0:
            for idx, align in enumerate(read_aligns):
                bam_discard.write(align)
                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 0, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))
            continue

        for idx, align in enumerate(read_aligns):
            if idx in keep:
                bam_keep.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 1, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))
            else:
                bam_discard.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 0, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))


def cluster_reads(input_bam: str, keep_bam: str, discard_bam: str, trim: int, mapping_quality_cutoff: int, alignment_stats: Optional[str] = None) -> Tuple[int, int, int, int]:


    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, 'wb', template=bam_in)
    bam_discard = AlignmentFile(discard_bam, 'wb', template=bam_in)

    if alignment_stats != None:
        alignment_stats_out = open(alignment_stats,'w')
        alignment_stats_out.write('read_id,mapping_id,filter_retained,query_start,query_end,mapq\n')


    num_reads = 0
    num_reads_kept = 0
    num_aligns = 0
    num_aligns_kept = 0
    for read_aligns in read_mappings_iter(bam_in):
        #read_aligns is a list of pysam.AlignmentFile objects
        num_reads += 1
        num_aligns += len(read_aligns)

        #non_substring_aligns = [aligns[x] for x in substring_filter(read_aligns,mapping_quality_cutoff = mapping_quality_cutoff)]

        keep = cluster_aligned_segments(read_aligns, trim, mapping_quality_cutoff)
        #keep is a list of indices into the read_aligns object of which alignments to keep

        if len(keep) == 0:
            for idx, align in enumerate(read_aligns):
                bam_discard.write(align)
                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 0, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))

            continue

        for idx, align in enumerate(read_aligns):
            if idx in keep:
                bam_keep.write(read_aligns[idx])

                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 1, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))
            else:
                bam_discard.write(read_aligns[idx])
                if alignment_stats != None:
                    alignment_stats_out.write('{read_id},{mapping_id},{filter_retained},{q_start},{q_end},{mapq}\n'.format(read_id = align.query_name, mapping_id = idx, filter_retained = 0, q_start = align.query_alignment_start, q_end = align.query_alignment_end, mapq = align.mapq))

        num_aligns_kept += len(keep)
    return (num_reads, num_reads_kept, num_aligns, num_aligns_kept)
