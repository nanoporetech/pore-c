import bisect
import gzip
import os
import sys
from collections import Counter, defaultdict, namedtuple
from dataclasses import dataclass
from enum import Enum, unique
from typing import (Dict, Iterator, List, NamedTuple, Optional, Pattern, Tuple,
                    Union)

import numpy as np
from pysam import AlignedSegment, AlignmentFile

from pore_c.model import FragmentMap


@unique
class AlignedSegmentMappingType(Enum):
    simple = "simple"
    multi_frag = "multi_frag"


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

    def assign_to_fragment(self, frag_map: "FragmentMap") -> None:
        """Find overlapping fragments and assign this segment to a single fragment.

        If an aligned segment overlaps multiple fragments then we take the longest overlap as
        the canonical one [TODO: should we take the one with highest MQ?]
        """
        ###TODO: this line needs to be replaced with something simpler, as the fragment assignments are already calculated:
        overlaps = list(frag_map.iter_overlaps((self.chrom, self.start, self.end)))
        num_overlaps = len(overlaps)
        if num_overlaps == 0:
            raise ValueError(
                "No overlap found, is input formatted correctly?. {}".format(self)
            )

        idx, mapping_type = None, None
        if num_overlaps == 1:
            # simple overlap
            mapping_type = AlignedSegmentMappingType.simple
        else:
            # case of multiple overlaps, take longest
            num_distinct_frags = len(set([o.frag_id for o in overlaps]))
            overlaps.sort(key=lambda x: x.overlap, reverse=True)
            if num_distinct_frags == 1:
                # shouldn't have multiple overlaps to a single fragment
                raise ValueError(self)
            else:
                # alignment crosses fragment boundary, pick longest overlap as fragment
                mapping_type = AlignedSegmentMappingType.multi_frag
        best_hit = overlaps[0]
        #        print('from overlaps: {}\nThis overlap was chosen: {}'.format(overlaps, best_hit))
        self.set_fragment(
            best_hit.frag_id, best_hit.frag_end, best_hit.overlap, mapping_type
        )

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
            align.mapping_quality,
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
            strand="16" if self.strand else "0",
            chrom=self.chrom,
            frag_endpoint=self.frag_endpoint,
            frag_id=self.frag_id,
        )


@dataclass
class ReadToFragments(object):
    """Holds the fragment assignments for a single read"""

    read_name: str
    num_frags: int
    num_nonadj_frags: int
    nonadj_vector: List[int]
    fragment_assignments: List[ReadToFragmentAssignment]

    def length(self):
        return len(self.fragment_assignments)

    def to_HiC_str(self):
        quals, mappings = zip(
            *[
                (str(_.max_mapping_quality), _.to_HiC_str())
                for _ in self.fragment_assignments
            ]
        )
        return "{name} {mappings} {quals}\n".format(
            name=self.read_name, mappings=" ".join(mappings), quals=" ".join(quals)
        )

    def stats(self):

        # percent of read mapped can be calculated once this table is combined with
        #   the fastq summary table, which contains the length of each read.
        tot_overlap = sum([_.total_overlap for _ in self.fragment_assignments])

        # intra:inter calculation
        intra = 0
        inter = 0

        for x1 in range(self.num_frags - 1):
            #skip a contact if the fragment wasn't adjacent to preceding fragment
            if not self.nonadj_vector[x1]:
                continue
            for x2 in range(x1 + 1, self.num_frags):
                #skip a contact if the fragment wasn't adjacent to preceding fragment
                if not self.nonadj_vector[x2]:
                    continue
                if (
                    self.fragment_assignments[x1].chrom
                    == self.fragment_assignments[x2].chrom
                ):
                    intra += 1
                else:
                    inter += 1


        uniqs = set([self.fragment_assignments[x].chrom if self.nonadj_vector[x] else "NonUniqueFragment" for x in range(self.num_frags)])

        uniqs.discard("NonUniqueFragment")


        if intra > 0 and inter > 0:
            intra_ratio = intra / float(inter)

        #all contacts are intra (inter == 0 is implied)
        elif intra > 0:
            intra_ratio = "Inf"
        elif inter > 0:
            intra_ratio = 0.0
        #no contacts at all
        else: 
            intra_ratio = "NaN"

        return "{name},{contact_count},{num_aligned_bases},{num_nonadj_frags},{uniq_chrs},{intra},{inter},{pct_intra},{intra_ratio}\n".format(
            name=self.read_name,
            contact_count=self.num_frags,
            num_aligned_bases=tot_overlap,
            num_nonadj_frags=self.num_nonadj_frags,
            uniq_chrs=len(uniqs),
            intra=intra,
            inter=inter,
            pct_intra=intra / float(intra + inter) if intra + inter > 0 else "NaN",
            intra_ratio=intra_ratio
        )

    def groupBy(self):
        self.overlaps = defaultdict(list)
        for align in self.fragment_assignments:
            self.overlaps[(align.read_name, align.read_start)].append(align)

    @classmethod
    def from_read_alignments(
        cls, read_aligns: "ReadAlignments", fragment_map: FragmentMap
    ):
        #fragment ordering in the read isn't preserved by dictionary
        frag_overlap_order = []
        frag_overlaps = defaultdict(list)
        # read_aligns.aligns is now a list of BedHits, which have been selected from overlapping reads already
        for idx, align in enumerate( read_aligns.aligns ):
            if align.frag_id not in frag_overlaps:
                frag_overlap_order.append(align.frag_id)
            frag_overlaps[align.frag_id].append(align)
          
        num_frags = len(frag_overlaps)
        print(frag_overlaps.keys(), frag_overlap_order)
        assert num_frags == len(frag_overlap_order)
        if num_frags > 1:
            #this list preserves the order of observation within the read.
            frag_ids = frag_overlap_order
#            frag_ids = list(frag_overlaps.keys())
            # FIXME: edge case where telomeric fragments from different chromosomes considered adjacent
            # n.b.: this is the number of non-adjacent pairs. in a concat a b c d e, i.e., there
            #      should be at most (a,b),(b,c),(c,d),(d,e): 4 pairs of potentially non-adjacent monomers
            #change this so that the nonadj vector is a boolean of whether a frag is to be included 
            # in uniqueness calcs such as cis:trans. First frag always included, everything subsequent
            # is identical to the old non-adjacency list that was formerly being summed.
            nonadj_vector = [True] + [(abs(b - a) > 1) for a, b in zip(frag_ids[:-1], frag_ids[1:])]
            num_nonadj_frags = sum(nonadj_vector) - 1
            ###the appropriate behaviour here would be if a pair of fragments are adjacent, one of the fragments should be removed
        else:
            num_nonadj_frags = 0
            nonadj_vector = [True]
        fragment_assignments = []
        for idx, data in enumerate( frag_overlaps.items()):
            frag_id, aligns = data
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
                sum(overlap_lengths),
            )
            if nonadj_vector[idx]:
                fragment_assignments.append(r)
        num_frags = sum(nonadj_vector)
        return cls(
            read_aligns.read_name, num_frags, num_nonadj_frags, nonadj_vector, fragment_assignments
        )


@dataclass
class BedHit(object):
    # format: {ch}\t{t_st}\t{t_en}\t{strand}\t{read_id}\t{q_st}\t{q_en}\t{mapq}\t{frag_ch}\t{frag_st}\t{frag_en}\t{frag_id}\t{overlap}
    chrom: str
    start: int
    end: int
    strand: bool
    read_name: str
    read_start: int
    read_end: int
    mapping_quality: int
    frag_ch: str
    frag_start: int
    frag_end: int
    frag_id: int
    frag_overlap: int

    @classmethod
    def from_bedformat(cls, align):
        l = align.strip().split()
        for field in [1, 2, 3, 5, 6, 7, 9, 10, 11, 12]:
            l[field] = int(l[field])
        return cls(*l)


@dataclass
class ReadAlignments(object):
    """Holds the aligned segments for a single read"""

    read_name: str
    aligns: List[AlignedSegmentToFragment]
    #    aligns: Dict #keyed on query_start, value of a list of BedHit objects, sorted by

    @staticmethod
    def iter_bam(input_bam: str) -> Iterator["ReadAlignments"]:
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
                    raise IOError(
                        "Seen this read already, is the BAM namesorted?: {}".format(
                            align.read_name
                        )
                    )
                yield ReadAlignments(
                    current_read_name, sorted(aligns, key=lambda x: x.read_start)
                )
                reads_seen.add(current_read_name)
                aligns = [align]
                current_read_name = align.read_name
        yield ReadAlignments(
            current_read_name, sorted(aligns, key=lambda x: x.read_start)
        )

    @staticmethod
    def iter_bed(input_bed: str) -> Iterator["ReadAlignments"]:
        """Iterate over a namesorted Bam and extract all of the aligments for a given read"""
        pre_aligns = defaultdict(list)
        aligns = []
        current_read_name = None
        reads_seen = set([])
        for align in map(BedHit.from_bedformat, open(input_bed)):
            #            print("read name:",align.read_name)
            # format: {ch}\t{t_st}\t{t_en}\t{strand}\t{read_id}\t{q_st}\t{q_en}\t{mapq}\t{frag_ch}\t{frag_st}\t{frag_en}\t{frag_id}\t{overlap}
            if current_read_name is None:
                current_read_name = align.read_name
                pre_aligns[align.read_start].append(align)
            elif current_read_name == align.read_name:
                pre_aligns[align.read_start].append(align)
            else:
                if align.read_name in reads_seen:
                    raise IOError(
                        "Seen this read already in set {}, is the BAM namesorted?: {}".format(
                            reads_seen, align.read_name
                        )
                    )
                aligns = []
                for overlaps in pre_aligns.values():
                    sorted_overlaps = sorted(overlaps, key=lambda x: x.frag_overlap)
                    aligns.append(sorted_overlaps[-1])
                yield ReadAlignments(
                    current_read_name, sorted(aligns, key=lambda x: x.read_start)
                )
                reads_seen.add(current_read_name)
                pre_aligns = defaultdict(list)
                pre_aligns[align.read_start].append(align)
                current_read_name = align.read_name
        yield ReadAlignments(
            current_read_name, sorted(aligns, key=lambda x: x.read_start)
        )


def map_to_fragments(
    input_bam: str,
    bed_file: str,
    output_file: str,
    method: str,
    stats_file: Optional[str] = None,
) -> None:
    fm = FragmentMap.from_bed_file(bed_file)
    f_out = open(output_file, "w")
    if stats_file != None:
        stats_out = open(stats_file, "w")
        stats_out.write(
            "read_id,contact_count,num_aligned_bases,num_nonadj_frags,uniq_chrs,intra,inter,pct_intra,intra_ratio\n"
        )
    for read_alignments in ReadAlignments.iter_bed(input_bam):
        frag_mapping = ReadToFragments.from_read_alignments(read_alignments, fm)
        # do not print out entries to a .poreC file for a single monomer. An entry must be at least a pairwise contact.
        if frag_mapping.length() > 1:
            f_out.write(frag_mapping.to_HiC_str())
        if stats_file != None:
            stats_out.write(frag_mapping.stats())

    f_out.close()
    if stats_file:
        stats_out.close()
