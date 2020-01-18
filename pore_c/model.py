from enum import Enum
from typing import Dict, List, NewType, Optional

import numpy as np
import pandas as pd
import pysam
from pydantic import BaseModel, confloat, conint, constr

from .config import (
    ALIGN_IDX_DTYPE,
    ALIGN_SCORE_DTYPE,
    FRAG_IDX_DTYPE,
    GENOMIC_COORD_DTYPE,
    GENOMIC_DISTANCE_DTYPE,
    HAPLOTYPE_IDX_DTYPE,
    MQ_DTYPE,
    PERCENTAGE_DTYPE,
    PHASE_SET_DTYPE,
    READ_COORD_DTYPE,
    READ_DISTANCE_DTYPE,
    READ_IDX_DTYPE,
    STRAND_DTYPE,
)
from .utils import mean_qscore


class _BaseModel(BaseModel):
    @classmethod
    def pandas_dtype(cls, overrides=None):
        res = {}
        if overrides is None:
            overrides = {}
        overrides = overrides if overrides is not None else {}
        for column, col_schema in cls.schema()["properties"].items():
            if column in overrides:
                res[column] = overrides[column]
            else:
                dtype = col_schema.get("dtype")
                if dtype == "category" and "enum" in col_schema:
                    dtype = pd.CategoricalDtype(col_schema["enum"], ordered=True)
                res[column] = dtype
        return res

    def to_tuple(self):
        return tuple([_[1] for _ in self])

    @classmethod
    def to_dataframe(cls, data: List, overrides=Optional[Dict]):
        columns = [a[0] for a in data[0]]
        dtype = cls.pandas_dtype(overrides=overrides)
        df = pd.DataFrame([a.to_tuple() for a in data], columns=columns).astype(dtype)
        return df


class AlignmentType(str, Enum):
    unmapped = "unmapped"
    primary = "primary"
    secondary = "secondary"
    supplementary = "supplementary"


class FragmentRecord(_BaseModel):
    """Meta-data associated with a restriction fragments"""

    chrom: constr(min_length=1, strip_whitespace=True)
    start: conint(ge=0)
    end: conint(ge=0)
    fragment_id: conint(ge=1, strict=True)
    fragment_length: conint(ge=1, strict=True)

    class Config:
        use_enum_values = True
        fields = dict(
            chrom=dict(description="The chromosome/contig the fragment is derived from", dtype="category"),
            start=dict(
                description="The zero-based start position on the genome of the fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            end=dict(
                description="The zero-based end position on the genome of the fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            fragment_id=dict(description="Unique integer ID of the fragment, starts at 1", dtype=FRAG_IDX_DTYPE),
            fragment_length=dict(description="Length of the fragment", dtype=GENOMIC_COORD_DTYPE),
        )


class AlignmentRecord(_BaseModel):
    """A subset of the fields in the BAM file"""

    read_idx: conint(ge=0, strict=True)
    align_idx: conint(ge=0, strict=True)
    align_type: AlignmentType
    chrom: constr(min_length=1, strip_whitespace=True)
    start: conint(ge=0)
    end: conint(ge=0)
    strand: STRAND_DTYPE
    read_name: constr(min_length=1, strip_whitespace=True)
    read_length: conint(ge=1)
    read_start: conint(ge=0)
    read_end: conint(ge=0)
    mapping_quality: conint(ge=0, le=255)
    align_score: conint(ge=0)
    align_base_qscore: conint(ge=0)
    phase_set: int = 0
    phase_qual: conint(ge=0) = 0
    haplotype: conint(ge=-1) = -1

    class Config:
        use_enum_values = True
        fields = dict(
            read_idx=dict(description="Unique integer ID of the read", dtype=READ_IDX_DTYPE),
            align_idx=dict(description="Unique integer ID of the aligned segment", dtype=ALIGN_IDX_DTYPE),
            align_type=dict(description="The type of alignment", dtype="category"),
            chrom=dict(description="The chromosome/contig the read is aligned to", dtype="category"),
            start=dict(
                description="The zero-based start position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE
            ),
            end=dict(description="The end position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE),
            strand=dict(description="The alignment strand", dtype="bool"),
            read_name=dict(description="The original read name", dtype="str"),
            read_length=dict(description="The length of the read in bases", dtype=READ_COORD_DTYPE),
            read_start=dict(description="The start coordinate on the read (0-based)", dtype=READ_COORD_DTYPE),
            read_end=dict(description="The end coordinate on the read (0-based)", dtype=READ_COORD_DTYPE),
            mapping_quality=dict(description="The mapping quality as calculated by the aligner", dtype=MQ_DTYPE),
            align_score=dict(description="The alignment score as calculated by the aligner", dtype=ALIGN_SCORE_DTYPE),
            align_base_qscore=dict(
                description="The mean read base score for the aligned segment (rounded to the nearest integer).",
                dtype=ALIGN_SCORE_DTYPE,
            ),
            phase_set=dict(
                description="The ID of the phase set, often this is the start position of the phase block",
                dtype=PHASE_SET_DTYPE,
            ),
            phase_qual=dict(description="The phred-scaled quality score of the haplotype assignment", dtype=MQ_DTYPE),
            haplotype=dict(
                description=(
                    "The id of the haplotype within this block, usually set to 1 or 2. "
                    "A value of -1 means that this alignment is unphased"
                ),
                dtype=HAPLOTYPE_IDX_DTYPE,
            ),
        )

    @classmethod
    def from_aligned_segment(cls, align: pysam.AlignedSegment) -> "AlignmentRecord":
        """Extract information from a pysam Aligned segment"""
        read_name, read_idx, align_idx = align.query_name.split(":")
        read_idx, align_idx = int(read_idx), int(align_idx)

        if align.is_unmapped:
            align_cat = "unmapped"
            chrom, start, end, align_score = "NULL", 0, 0, 0
            read_length = align.query_length
            quals = align.query_qualities
            # TODO: handle this more gracefully
            if quals is None:
                align_base_qscore = 0
            else:
                align_base_qscore = mean_qscore(np.array(align.query_qualities))
        else:
            chrom, start, end = (align.reference_name, align.reference_start, align.reference_end)
            read_length = align.infer_read_length()
            align_score = align.get_tag("AS")
            align_base_qscore = mean_qscore(np.array(align.query_alignment_qualities))
            if align.is_secondary:
                align_cat = "secondary"
            elif align.is_supplementary:
                align_cat = "supplementary"
            else:
                align_cat = "primary"

        optional = {}
        for key, tag in [("haplotype", "HP"), ("phase_set", "PS"), ("phase_qual", "PC")]:
            if align.has_tag(tag):
                optional[key] = int(align.get_tag(tag))
        return cls(
            read_idx=read_idx,
            align_idx=align_idx,
            align_type=align_cat,
            chrom=chrom,
            start=start,
            end=end,
            strand=not align.is_reverse,
            read_name=read_name,
            read_length=read_length,
            read_start=align.query_alignment_start,
            read_end=align.query_alignment_end,
            mapping_quality=align.mapq,
            align_score=align_score,
            align_base_qscore=np.rint(align_base_qscore),
            **optional
        )

    @classmethod
    def to_dataframe(cls, aligns: List, chrom_order: List[str] = None):
        columns = [a[0] for a in aligns[0]]
        if chrom_order:
            overrides = {"chrom": pd.CategoricalDtype(chrom_order, ordered=True)}
        else:
            overrides = {}
        dtype = cls.pandas_dtype(overrides=overrides)
        df = pd.DataFrame([a.to_tuple() for a in aligns], columns=columns).astype(dtype)
        return df


class AlignmentFilterReason(str, Enum):
    null = "null"
    Pass = "pass"
    unmapped = "unmapped"
    singleton = "singleton"
    low_mq = "low_mq"
    overlap_on_read = "overlap_on_read"
    not_on_shortest_path = "not_on_shortest_path"


class PoreCRecord(AlignmentRecord):
    pass_filter: bool = True
    filter_reason: AlignmentFilterReason = AlignmentFilterReason.null
    fragment_id: conint(ge=0) = 0
    num_contained_fragments: conint(ge=0) = 0
    num_overlapping_fragments: conint(ge=0) = 0
    overlap_length: conint(ge=0) = 0
    fragment_start: conint(ge=0) = 0
    fragment_end: conint(ge=0) = 0
    perc_of_alignment: confloat(ge=0, le=100) = 0.0
    perc_of_fragment: confloat(ge=0, le=100) = 0.0
    is_contained: bool = False

    class Config:
        use_enum_values = True
        fields = dict(
            pass_filter=dict(description="Boolean flag, true if alignment passes all filters", dtype="bool"),
            filter_reason=dict(
                description="If an alignment fails the filter the reason will be listed here", dtype="category"
            ),
            fragment_id=dict(
                description="The UID of the restriction fragment assigned to this alignment", dtype=FRAG_IDX_DTYPE
            ),
            num_contained_fragments=dict(
                description="The number of restriction fragments completely contained within this alignment",
                dtype="uint32",
            ),
            num_overlapping_fragments=dict(
                description="The number of restriction fragments overlapping this alignment", dtype="uint32"
            ),
            overlap_length=dict(
                description="The length of the overlap between alignment and fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            fragment_start=dict(
                description="The start point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            fragment_end=dict(
                description="The end point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            perc_of_alignment=dict(
                description="The percentage of the aligned segment that overlaps the assigned fragment",
                dtype=PERCENTAGE_DTYPE,
            ),
            perc_of_fragment=dict(
                description="The percentage of the assigned restriction fragment that overlaps the aligned segment",
                dtype=PERCENTAGE_DTYPE,
            ),
            is_contained=dict(
                description="Boolean flag to inidicate if the alignment is fully contained with the fragment",
                dtype="bool",
            ),
        )

    @classmethod
    def init_dataframe(cls, align_df: "AlignmentRecordDf") -> "PoreCRecordDf":
        res = align_df.copy()
        schema = cls.schema()["properties"]
        dtype = cls.pandas_dtype()
        additional_fields = set(dtype.keys()) - (set(align_df.index.names) | set(align_df.columns))
        num_rows = len(res)
        for column in [c for c in schema if c in additional_fields]:
            default_value = schema[column]["default"]
            res[column] = pd.Series([default_value] * num_rows, index=res.index).astype(dtype[column])
        return res


class HaplotypePairType(str, Enum):
    null = "null"
    trans = "trans"
    unphased = "unphased"
    semi_phased = "semi_phased"
    phased_sets_differ = "phased_sets_differ"
    phased_h_cis = "phased_h_cis"
    phased_h_trans = "phased_h_trans"


class PoreCContactRecord(_BaseModel):
    read_idx: conint(ge=0, strict=True)
    contact_is_direct: bool = False
    contact_is_cis: bool = False
    contact_read_distance: int = 0
    contact_genome_distance: int = 0
    contact_fragment_distance: conint(ge=0, strict=True)
    haplotype_pair_type: HaplotypePairType = HaplotypePairType.null
    align1_align_idx: conint(ge=0, strict=True)
    align1_chrom: constr(min_length=1, strip_whitespace=True)
    align1_start: conint(ge=0)
    align1_end: conint(ge=0)
    align1_strand: STRAND_DTYPE
    align1_mapping_quality: conint(ge=0, le=255)
    align1_align_score: conint(ge=0)
    align1_align_base_qscore: conint(ge=0)
    align1_phase_set: int = 0
    align1_phase_qual: int = 0
    align1_haplotype: conint(ge=-1) = -1
    align1_fragment_id: conint(ge=0) = 0
    align1_fragment_start: conint(ge=0) = 0
    align1_fragment_end: conint(ge=0) = 0
    align2_align_idx: conint(ge=0, strict=True)
    align2_chrom: constr(min_length=1, strip_whitespace=True)
    align2_start: conint(ge=0)
    align2_end: conint(ge=0)
    align2_strand: STRAND_DTYPE
    align2_mapping_quality: conint(ge=0, le=255)
    align2_align_score: conint(ge=0)
    align2_align_base_qscore: conint(ge=0)
    align2_phase_set: int = 0
    align1_phase_qual: int = 0
    align2_haplotype: conint(ge=-1) = -1
    align2_fragment_id: conint(ge=0) = 0
    align2_fragment_start: conint(ge=0) = 0
    align2_fragment_end: conint(ge=0) = 0

    class Config:
        use_enum_values = True
        fields = dict(
            read_idx=dict(description="Unique integer ID of the read", dtype=READ_IDX_DTYPE),
            contact_is_direct=dict(
                description="There are no intervening assigned restriction fragments on the read", dtype="bool"
            ),
            contact_is_cis=dict(description="Both alignments come from the same chromsome/contig", dtype="bool"),
            contact_read_distance=dict(
                description=(
                    "The distance between the end of the left alignment and the start of the right "
                    "alignment on the read"
                ),
                dtype=READ_DISTANCE_DTYPE,
            ),
            contact_genome_distance=dict(
                description=(
                    "The distance between the end of the left alignment and the start of the right alignment "
                    "(valid for cis contacts only)"
                ),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
            contact_fragment_distance=dict(
                description=(
                    "The distance between the midpoints of the assigned fragments (valid for cis contacts only)"
                ),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
            haplotype_pair_type=dict(
                description=(
                    "A categorical variable describing the relationship between the haplotypes assigned to each of the "
                    "alignments in a contact",
                ),
                dtype="category",
            ),
            align1_align_idx=dict(description="Unique integer ID of the first aligned segment", dtype=ALIGN_IDX_DTYPE),
            align1_chrom=dict(description="The chromosome/contig of the first aligned segment", dtype="category"),
            align1_start=dict(
                description="The zero-based start position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE
            ),
            align1_end=dict(description="The end position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE),
            align1_strand=dict(description="The alignment strand", dtype="bool"),
            align1_mapping_quality=dict(description="The mapping quality as calculated by the aligner", dtype=MQ_DTYPE),
            align1_align_score=dict(
                description="The alignment score as calculated by the aligner", dtype=ALIGN_SCORE_DTYPE
            ),
            align1_align_base_qscore=dict(
                description="The mean read base score for the aligned segment (rounded to the nearest integer).",
                dtype=ALIGN_SCORE_DTYPE,
            ),
            align1_phase_set=dict(
                description="The ID of the phase set, often this is the start position of the phase block",
                dtype=PHASE_SET_DTYPE,
            ),
            align1_phase_qual=dict(
                description="The phred-scaled quality score of the haplotype assignment", dtype=MQ_DTYPE
            ),
            align1_haplotype=dict(
                description=(
                    "The id of the haplotype within this block, usually set to 1 or 2. "
                    "A value of -1 means that this alignment is unphased"
                ),
                dtype=HAPLOTYPE_IDX_DTYPE,
            ),
            align1_fragment_id=dict(
                description="The UID of the restriction fragment assigned to this alignment", dtype=FRAG_IDX_DTYPE
            ),
            align1_fragment_start=dict(
                description="The start point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            align1_fragment_end=dict(
                description="The end point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            align2_align_idx=dict(description="Unique integer ID of the first aligned segment", dtype=ALIGN_IDX_DTYPE),
            align2_chrom=dict(description="The chromosome/contig of the first aligned segment", dtype="category"),
            align2_start=dict(
                description="The zero-based start position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE
            ),
            align2_end=dict(description="The end position on the genome of the alignment", dtype=GENOMIC_COORD_DTYPE),
            align2_strand=dict(description="The alignment strand", dtype="bool"),
            align2_mapping_quality=dict(description="The mapping quality as calculated by the aligner", dtype=MQ_DTYPE),
            align2_align_score=dict(
                description="The alignment score as calculated by the aligner", dtype=ALIGN_SCORE_DTYPE
            ),
            align2_align_base_qscore=dict(
                description="The mean read base score for the aligned segment (rounded to the nearest integer).",
                dtype=ALIGN_SCORE_DTYPE,
            ),
            align2_phase_set=dict(
                description="The ID of the phase set, often this is the start position of the phase block",
                dtype=PHASE_SET_DTYPE,
            ),
            align2_phase_qual=dict(
                description="The phred-scaled quality score of the haplotype assignment", dtype=MQ_DTYPE
            ),
            align2_haplotype=dict(
                description=(
                    "The id of the haplotype within this block, usually set to 1 or 2. "
                    "A value of -1 means that this alignment is unphased"
                ),
                dtype=HAPLOTYPE_IDX_DTYPE,
            ),
            align2_fragment_id=dict(
                description="The UID of the restriction fragment assigned to this alignment", dtype=FRAG_IDX_DTYPE
            ),
            align2_fragment_start=dict(
                description="The start point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
            align2_fragment_end=dict(
                description="The end point on the genome of this restriction fragment", dtype=GENOMIC_COORD_DTYPE
            ),
        )

    @classmethod
    def from_pore_c_align_pair(cls, read_idx: int, align1, align2, contact_is_direct: bool = False):
        contact_read_distance = align2.read_start - align1.read_end
        if align1.fragment_id > align2.fragment_id:
            align1, align2 = align2, align1

        contact_is_cis = align1.chrom == align2.chrom
        if contact_is_cis:
            contact_genome_distance = align2.start - align1.end
            contact_fragment_distance = align2.fragment_midpoint - align1.fragment_midpoint
        else:
            contact_genome_distance = 0
            contact_fragment_distance = 0

        haplotype_pair_type = HaplotypePairType.null
        if not contact_is_cis:
            haplotype_pair_type = HaplotypePairType.trans
        elif align1.haplotype == -1 and align2.haplotype == -1:
            haplotype_pair_type = HaplotypePairType.unphased
        elif align1.haplotype == -1 or align2.haplotype == -1:
            haplotype_pair_type = HaplotypePairType.semi_phased
        elif align1.phase_set != align2.phase_set:
            haplotype_pair_type = HaplotypePairType.phased_sets_differ
        elif align1.haplotype == align2.haplotype:
            haplotype_pair_type = HaplotypePairType.phased_h_cis
        else:
            haplotype_pair_type = HaplotypePairType.phased_h_trans

        return cls(
            read_idx=read_idx,
            contact_is_direct=contact_is_direct,
            contact_is_cis=contact_is_cis,
            contact_read_distance=contact_read_distance,
            contact_genome_distance=contact_genome_distance,
            contact_fragment_distance=contact_fragment_distance,
            haplotype_pair_type=haplotype_pair_type,
            align1_align_idx=align1.align_idx,
            align1_chrom=align1.chrom,
            align1_start=align1.start,
            align1_end=align1.end,
            align1_strand=align1.strand,
            align1_mapping_quality=align1.mapping_quality,
            align1_align_score=align1.align_score,
            align1_align_base_qscore=align1.align_base_qscore,
            align1_phase_set=align1.phase_set,
            align1_phase_qual=align1.phase_qual,
            align1_haplotype=align1.haplotype,
            align1_fragment_id=align1.fragment_id,
            align1_fragment_start=align1.fragment_start,
            align1_fragment_end=align1.fragment_end,
            align2_align_idx=align2.align_idx,
            align2_chrom=align2.chrom,
            align2_start=align2.start,
            align2_end=align2.end,
            align2_strand=align2.strand,
            align2_mapping_quality=align2.mapping_quality,
            align2_align_score=align2.align_score,
            align2_align_base_qscore=align2.align_base_qscore,
            align2_phase_set=align2.phase_set,
            align2_phase_qual=align2.phase_qual,
            align2_haplotype=align2.haplotype,
            align2_fragment_id=align2.fragment_id,
            align2_fragment_start=align2.fragment_start,
            align2_fragment_end=align2.fragment_end,
        )


class PoreCConcatemerRecord(_BaseModel):
    read_idx: conint(ge=0, strict=True)
    indirect_contacts: conint(ge=0, strict=True)
    direct_contacts: conint(ge=0, strict=True)
    total_contacts: conint(ge=0, strict=True)
    read_order: conint(ge=0, strict=True)
    total_cis_contacts: conint(ge=0, strict=True)
    haplotype_phased_h_cis: conint(ge=0, strict=True)
    haplotype_phased_h_trans: conint(ge=0, strict=True)
    haplotype_phased_sets_differ: conint(ge=0, strict=True)
    haplotype_semi_phased: conint(ge=0, strict=True)
    haplotype_unphased: conint(ge=0, strict=True)
    max_indirect_contact_genome_distance: conint(ge=0, strict=True)
    max_direct_contact_genome_distance: conint(ge=0, strict=True)
    max_indirect_contact_fragment_distance: conint(ge=0, strict=True)
    max_direct_contact_fragment_distance: conint(ge=0, strict=True)

    class Config:
        use_enum_values = True
        fields = dict(
            read_idx=dict(description="Unique integer ID of the read", dtype=READ_IDX_DTYPE),
            read_order=dict(description="The number of monomers for this read", dtype="uint32"),
            total_contacts=dict(description="The total number of contacts for this read", dtype="uint32"),
            direct_contacts=dict(
                description="The total number direct (adjacent on read) contacts for this read", dtype="uint32"
            ),
            indirect_contacts=dict(
                description="The total number indirect (non-adjacent on read) contacts for this read", dtype="uint32"
            ),
            total_cis_contacts=dict(
                description="The total number of cis-contacts (direct + indirect) for this read", dtype="uint32"
            ),
            haplotype_unphased=dict(
                description="The number of cis contacts where both members of the pair are unphased", dtype="uint32"
            ),
            haplotype_semi_phased=dict(
                description="The number of cis contacts where one member of the pair is unphased", dtype="uint32"
            ),
            haplotype_phased_sets_differ=dict(
                description=(
                    "The number of cis contacts where both members of the pair are phased but the phase sets differ"
                ),
                dtype="uint32",
            ),
            haplotype_phased_h_trans=dict(
                description=(
                    "The number of cis contacts where both members of the pair are phased, are part of the same phase "
                    "group, but the haplotypes differ"
                ),
                dtype="uint32",
            ),
            haplotype_phased_h_cis=dict(
                description=(
                    "The number of cis contacts where both members of the pair are phased, are part of the same phase "
                    "group, and the haplotypes agree"
                ),
                dtype="uint32",
            ),
            max_direct_contact_fragment_distance=dict(
                description=("The longest distance between fragment midpoints for all direct contacts",),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
            max_indirect_contact_fragment_distance=dict(
                description=("The longest distance between fragment midpoints for all indirect contacts",),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
            max_direct_contact_genome_distance=dict(
                description=("The longest distance between alignment endpoints for all direct contacts",),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
            max_indirect_contact_genome_distance=dict(
                description=("The longest distance between alignment endpoints for all indirect contacts",),
                dtype=GENOMIC_DISTANCE_DTYPE,
            ),
        )


AlignmentRecordDf = NewType("AlignmentRecordDf", pd.DataFrame)
FragmentRecordDf = NewType("FragmentRecordDf", pd.DataFrame)
PoreCRecordDf = NewType("PoreCRecordDf", pd.DataFrame)
PoreCContactRecordDf = NewType("PoreCContactRecordDf", pd.DataFrame)
PoreCConcatemerRecordDf = NewType("PoreCConcatemerRecordDf", pd.DataFrame)

Chrom = NewType("Chrom", str)
