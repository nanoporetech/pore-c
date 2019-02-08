from typing import Pattern, List, NamedTuple, Iterator, Union, Dict, Tuple
from pysam import AlignmentFile, AlignedSegment
import bisect
import sys
import gzip
import numpy as np
from collections import defaultdict

from enum import Enum, unique


from collections import namedtuple, Counter
from dataclasses import dataclass

from pore_c.model import FragmentMap


@unique
class AlignedSegmentMappingType(Enum):
    simple = 'simple'
    multi_frag = 'multi_frag'


@dataclass
class AlignedSegmentToFragment(object):
    """Represents a subset of the data in a pysam.AlignedSegment along with an assignment to a fragment"""
    chrom: str
    start: int
    end: int
    strand: bool
    read_name: str
    read_start: int
    read_end: int
    mapping_quality: int

    frag_id: int = -1
    frag_end: int = -1
    frag_overlap: int = -1
    mapping_type: AlignedSegmentMappingType = None

    def set_fragment(self, frag_id, frag_end, frag_overlap, mapping_type):
        self.frag_id = frag_id
        self.frag_end = frag_end
        self.frag_overlap = frag_overlap
        self.mapping_type = mapping_type

    def assign_to_fragment(self, frag_map: 'FragmentMap') -> None:
        """Find overlapping fragments and assign this segment to a single fragment.

        If an aligned segment overlaps multiple fragments then we take the longest overlap as
        the canonical one [TODO: should we take the one with highest MQ?]
        """
        overlaps = list(frag_map.iter_overlaps((self.chrom, self.start, self.end)))
        num_overlaps = len(overlaps)
        if num_overlaps == 0:
            raise ValueError("No overlap found, is input formatted correctly?. {}".format(self))

        idx, mapping_type = None, None
        if num_overlaps == 1:
            #simple overlap
            mapping_type = AlignedSegmentMappingType.simple
        else:
            # case of multiple overlaps, take longest
            num_distinct_frags = len(set([o.frag_id for o in overlaps]))
            overlaps.sort(key=lambda x: x.overlap, reverse=True)
            if num_distinct_frags == 1:
                # shouldn't have multiple overlaps to a single fragment
                raise ValueError(self)
            else:
                #alignment crosses fragment boundary, pick longest overlap as fragment
                mapping_type = AlignedSegmentMappingType.multi_frag
        best_hit = overlaps[0]
        self.set_fragment(best_hit.frag_id, best_hit.frag_end, best_hit.overlap, mapping_type)


    @classmethod
    def from_pysam(cls, align):
        return cls(
            align.reference_name,
            align.reference_start,
            align.reference_end,
            align.is_reverse,
            align.query_name,
            align.query_alignment_start,
            align.query_alignment_end,
            align.mapping_quality
        )


@dataclass
class ReadToFragmentAssignment(object):
    """Holds data about how well a read maps to a single fragment"""
    frag_id: int
    frag_endpoint: int
    chrom: str
    strand: bool
    num_aligns: int
    max_mapping_quality: int
    best_overlap: int
    total_overlap: int

    def to_HiC_str(self):
        return "{strand} {chrom} {frag_endpoint} {frag_id}".format(
            strand = '16' if self.strand else 0,
            chrom = self.chrom, frag_endpoint = self.frag_endpoint, frag_id=self.frag_id
        )


@dataclass
class ReadToFragments(object):
    """Holds the fragment assignments for a single read"""
    read_name: str
    num_frags: int
    num_nonadj_frags: int
    fragment_assignments: List[ReadToFragmentAssignment]

    def to_HiC_str(self):
        quals, mappings = zip(*[
            (str(_.max_mapping_quality), _.to_HiC_str())
            for _ in self.fragment_assignments
        ])
        return "{name} {mappings} {quals}\n".format(
            name=self.read_name, mappings = " ".join(mappings), quals= " ".join(quals)
        )

    def log(self):
        #percent of read mapped can be calculated once this table is combined with
        #   the fastq summary table, which contains the length of each read.
        tot_overlap = sum([_.overlap for _ in self.fragment_assignments])

        return "{name},{contact_count},{num_aligned_bases},{num_nonadj_frags}".format(
            name = self.read_name, contact_count = self.num_frags,
            num_aligned_bases = tot_overlap, num_nonadj_frags = self.num_nonadj_frags)


    @classmethod
    def from_read_alignments(cls, read_aligns: 'ReadAlignments', fragment_map: FragmentMap):
        frag_overlaps = defaultdict(list)
        for idx, align in enumerate(read_aligns.aligns):
            align.assign_to_fragment(fragment_map)
            frag_overlaps[align.frag_id].append(align)
        num_frags = len(frag_overlaps)
        if num_frags > 1:
            frag_ids = sorted(frag_overlaps.keys())
            #FIXME: edge case where telomeric fragments from different chromosomes considered adjacent
            #n.b.: this is the number of non-adjacent pairs. in a concat a b c d e, i.e., there
            #      should be at most (a,b),(b,c),(c,d),(d,e): 4 pairs of potentially non-adjacent monomers
            num_nonadj_frags = sum(([(b - a > 1) for a, b in zip(frag_ids[:-1], frag_ids[1:])]))
        else:
            num_nonadj_frags = 0
        fragment_assignments = []
        for frag_id, aligns in frag_overlaps.items():
            aligns.sort(key=lambda x: x.mapping_quality, reverse=True)
            best_hit = aligns[0]
            overlap_lengths = [_.frag_overlap for _ in aligns]
            r = ReadToFragmentAssignment(
                frag_id,
                best_hit.frag_end,
                best_hit.chrom,
                best_hit.strand,
                len(aligns),
                best_hit.mapping_quality,
                best_hit.frag_overlap,
                sum(overlap_lengths)
            )
            fragment_assignments.append(r)
        return cls(read_aligns.read_name, num_frags, num_nonadj_frags, fragment_assignments)


@dataclass
class ReadAlignments(object):
    """Holds the aligned segments for a single read"""
    read_name: str
    aligns: List[AlignedSegmentToFragment]

    @staticmethod
    def iter_bam(input_bam: str) -> Iterator['ReadAlignments']:
        """Iterate over a namesorted Bam and extract all of the aligments for a given read"""
        aligns = []
        current_read_name = None
        reads_seen = set([])
        for align in map(AlignedSegmentToFragment.from_pysam, AlignmentFile(input_bam)):
            if current_read_name is None:
                current_read_name = align.read_name
                aligns.append(align)
            elif current_read_name == align.read_name:
                aligns.append(align)
            else:
                if align.read_name in reads_seen:
                    raise IOError("Seen this read already, is the BAM namesorted?: {}".format(align.read_name))
                yield ReadAlignments(current_read_name, sorted(aligns, key=lambda x: x.read_start))
                reads_seen.add(current_read_name)
                aligns = [align]
                current_read_name = align.read_name
        yield ReadAlignments(current_read_name, sorted(aligns, key=lambda x: x.read_start))


def map_to_fragments(input_bam: str, bed_file: str, output_file: str, method: str) -> None:
    fm = FragmentMap.from_bed_file(bed_file)
    f_out = open(output_file, 'w')
    for read_alignments in ReadAlignments.iter_bam(input_bam):
        frag_mapping = ReadToFragments.from_read_alignments(read_alignments, fm)
        f_out.write(frag_mapping.to_HiC_str())
    f_out.close()

def map_to_fragments_logging(input_bam: str, bed_file: str, output_file: str, method: str, log_file: str) -> None:
    fm = FragmentMap.from_bed_file(bed_file)
    f_out = open(output_file, 'w')
    log_out = open(log_file, 'w')
    log_out.write("name,contact_count,num_aligned_bases,num_nonadj_frags")
    for read_alignments in ReadAlignments.iter_bam(input_bam):
        frag_mapping = ReadToFragments.from_read_alignments(read_alignments, fm)
        f_out.write(frag_mapping.to_HiC_str())
        log_out.write(frag_mapping.log())

    f_out.close()
    log_out.close()
