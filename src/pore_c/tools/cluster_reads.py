import sys
from typing import Tuple, Optional
import numpy as np
from pysam import AlignmentFile
from intervaltree import IntervalTree

import random

def read_mappings_iter(bam):
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

def measure_overlaps(input_bam: str, output_table: str, no_zero):

    bam_in = AlignmentFile(input_bam)
    f_out = open(output_table,'w')
    header = "read_id,align1,align2,overlap_length\n"
    f_out.write(header)

    for read_aligns in read_mappings_iter(bam_in):
        l = len(read_aligns)
        for x in range(l-1):
            for y in range(1,l):
                olap = overlap(read_aligns[x],read_aligns[y])
                if no_zero and olap == 0:
                    continue
                f_out.write('{read_id},{align1},{align2},{overlap}\n'.format(read_id = read_aligns[x].query_name,align1 = x,align2 = y, overlap = olap))


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


##################
###fragdag code###
##################
def minimap_gapscore(length, o1=4, o2=24, e1=2, e2=1):
    return min([int(o1) + int(length) * int(e1), int(o2) + int(length) * int(e2)])

#O=5 E=2 is default for bwa bwasw
#O=6 E=1 is default for bwa mem

def bwa_gapscore(length, O = 5,E = 2):
    return (O + length * E)

#aligns must be sorted by END position (I think)
def fragDAG(_aligns, mapping_quality_cutoff = 0, aligner = "minimap2", params = "default"):
    #filter by mapq and (start,stop) first

    seen = {}
    for idx, align in enumerate(_aligns):
        if align.mapq >= mapping_quality_cutoff:
            coords = (align.query_alignment_start,align.query_alignment_emd) 
            if coords not in seen:
                seen[coords] = align.get_tag("AS")
            elif align.get_tag("AS") > seen[coords]:
                seen[coords] = align.get_tag("AS")
    keep = seen.keys()

    ##sort by end position so that DAG is preserved from cycles
    aligns = sorted([_aligns[x] for x in keep], key = lambda x: x.query_alignment_end)

    G = nx.DiGraph()
    edge_values = {}
    G.add_node("IN")
    G.add_node("OUT")
    L = len(aligns)
    readlen = aligns[0].query_length
    for x in keep:
        G.add_node(x)
        if aligner == "minimap2":
            GP_in = minimap_gapscore(aligns[x].query_alignment_start)
            GP_out = minimap_gapscore(readlen - aligns[x].query_alignment_end)
        elif aligner = "bwasw":
            GP_in = bwa_gapscore(aligns[x].query_alignment_start)
            GP_out = bwap_gapscore(readlen - aligns[x].query_alignment_end)
        AS = aligns[x].get_tag("AS")
        edge_values[("IN",x)] = (GP_in,AS) #store this data for diagnosis
        edge_values[(x,"OUT")] = (GP_out,AS) #store this data for diagnosis
        G.add_edge("IN",x,weight = GP_in - AS)
        G.add_edge(x,"OUT",weight = GP_out)
    for x in range(L-1):
        for y in range(x + 1, L):
            #if the two reads end on the same base, check which one starts first. If they map to the same interval, skip.
            # . otherwise, edge from the alignment that starts closer to the beginning of the read
            if aligns[x].query_alignment_end == aligns[y].query_alignment_end:
                if aligns[x].query_alignment_start == aligns[y].query_alignment_start:
                    continue 
                elif aligns[x].query_alignment_start > aligns[y].query_alignment_start:
                    align1 =  aligns[y]
                    align2 = aligns[x]
                else:
                    align1 =  aligns[x]
                    align2 = aligns[y]

            #abs, because if the alignments are overlapping, then this still describes the penalty of having to disregard the part
            # . of the second alignment that overlaps the first alignment
            if aligner == "minimap2":
                GP = minimap_gapscore(abs(align1.query_alignment_start - align2.query_alignment_end ))
            elif aligner = "bwasw":
                GP = bwa_gapscore(abs(aligns1.query_alignment_start - align2.query_alignment_end ))
            AS = align2.get_tag("AS")
            edge_values[(x,y)] = (GP,AS) #store this data for diagnosis
            G.add_edge(x,y, weight = GP - AS)

    #returns a list of the path through the graph as well as the 
    return nx.single_source_bellman_ford(G,"in", "out")[1]# edge_values # including this as output enables reconstruction of the DAG that ended up filtering reads out. useful for diagnostic purposes



def fragment_DAG(input_bam: str, keep_bam: str, discard_bam: str, mapping_quality_cutoff: int, aligner: str, aligner_params: str,  alignment_stats: Optional[str] = None):

    bam_in = AlignmentFile(input_bam)
    bam_keep = AlignmentFile(keep_bam, 'wb', template=bam_in)
    bam_discard = AlignmentFile(discard_bam, 'wb', template=bam_in)

    if alignment_stats != None:
        alignment_stats_out = open(alignment_stats,'w')
        alignment_stats_out.write('read_id,mapping_id,filter_retained,query_start,query_end,mapq\n')


    for read_aligns in read_mappings_iter(bam_in):

        #this line returns two lists, the number of the reads that are part of the best scoring path, and the scores for that path based on the scoring function outlined
        keep = fragDAG(read_aligns,mapping_quality_cutoff = mapping_quality_cutoff, aligner = aligner, params = aligner_params)

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
