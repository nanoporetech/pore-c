from enum import Enum
from typing import Dict, List, NewType, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
from ncls import NCLS
from pydantic import BaseModel, confloat, conint, constr
from pysam import AlignedSegment

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
                dtype='bool'
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
                description="There are no intervening assigned restriction fragments on the read", dtype='bool'
            ),
            contact_is_cis=dict(description="Both alignments come from the same chromsome/contig", dtype='bool'),
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


# GENOMIC_COORD_DTYPE = np.uint32  # should be fine as long as individual chromosomes are less than 4Gb
# READ_COORD_DTYPE = np.uint32
# FRAG_IDX_DTYPE = np.uint32
# READ_IDX_DTYPE = np.uint32
# ALIGN_IDX_DTYPE = np.uint32
# PERCENTAGE_DTYPE = np.float32
# HAPLOTYPE_IDX_DTYPE = np.int8


class basePorecDf(object):
    DTYPE = {}
    NANS = {}

    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        self._validation_errors = {}
        self._additional_columns = []
        self._check(pandas_obj)

    def _check(self, obj):
        # FIXFIX: this is pretty broken, currently need to cast before check
        # https://stackoverflow.com/questions/22697773/how-to-check-the-dtype-of-a-column-in-python-pandas
        errors = {}
        for key, value in self.DTYPE.items():
            _dtype = obj.dtypes.get(key, None)
            try:
                if _dtype is None:
                    errors[key] = "Missing column"
                elif value is None:
                    errors[key] = "Unset datatype"
                elif type(_dtype) is not type(value) and not isinstance(_dtype, value) and not _dtype == value:
                    errors[key] = "Mismatched dtype ({}/{}) (expected/found)".format(_dtype, value)
            except Exception as exp:
                raise ValueError("{}: {} - {}\n{}".format(key, value, _dtype, exp))
        self._validation_errors = errors
        for column, dtype in obj.dtypes.items():
            if column not in self.DTYPE:
                self._additional_columns.append(column)

    def assert_valid(self):
        if self._validation_errors:
            raise ValueError(
                "Failed validation:\n{}\n".format(
                    "\n".join(["{}: {}".format(*_) for _ in self._validation_errors.items()])
                )
            )

    @classmethod
    def set_dtype(cls, key, value):
        """If we don't know a datatype until execution time we can set it here"""
        cls.DTYPE[key] = value

    def cast(self, subset=False, fillna=False) -> pd.DataFrame:
        cols = list(self.DTYPE.keys())
        obj_cols = list(self._obj.columns)
        different_cols = set(cols).symmetric_difference(set(obj_cols))
        if different_cols:
            missing_cols = set(cols) - set(obj_cols)
            extra_cols = set(obj_cols) - set(cols)
            if missing_cols:
                raise ValueError("Columns missing from dataframe: {}".format(",".join(missing_cols)))
            if extra_cols:
                raise ValueError("Columns extrac in dataframe: {}".format(",".join(extra_cols)))
        type_df = {key: val for key, val in self.DTYPE.items() if val is not None}
        fillna_dict = self.NANS if fillna else {}

        _df = self._obj.loc[:, cols]
        _df = _df.fillna(fillna_dict)
        _df = _df.astype(type_df)
        return _df


@pd.api.extensions.register_dataframe_accessor("bamdf")
class BamEntryDf(basePorecDf):
    """An extension to handle a dataframe derived from a BAM file"""

    _ACCESSOR = "bamdf"
    DTYPE = {
        "read_idx": READ_IDX_DTYPE,
        "align_idx": ALIGN_IDX_DTYPE,
        "mapping_type": pd.CategoricalDtype(["unmapped", "primary", "supplementary", "secondary"], ordered=True),
        "chrom": str,
        "start": GENOMIC_COORD_DTYPE,
        "end": GENOMIC_COORD_DTYPE,
        "strand": bool,
        "read_name": str,
        "read_length": READ_COORD_DTYPE,
        "read_start": READ_COORD_DTYPE,
        "read_end": READ_COORD_DTYPE,
        "mapping_quality": np.uint8,
        "score": np.uint32,
    }

    @staticmethod
    def align_to_tuple(align: AlignedSegment) -> Tuple:
        read_name, read_idx, align_idx = align.query_name.split(":")
        read_idx, align_idx = int(read_idx), int(align_idx)
        if align.is_unmapped:
            align_cat = "unmapped"
            chrom, start, end, align_score = "NULL", 0, 0, 0
            read_length = align.query_length
        else:
            chrom, start, end = (align.reference_name, align.reference_start, align.reference_end)
            read_length = align.infer_read_length()
            align_score = align.get_tag("AS")
            if align.is_secondary:
                align_cat = "secondary"
            elif align.is_supplementary:
                align_cat = "supplementary"
            else:
                align_cat = "primary"
        return (
            read_idx,
            align_idx,
            align_cat,
            chrom,
            start,
            end,
            not align.is_reverse,
            read_name,
            read_length,
            align.query_alignment_start,
            align.query_alignment_end,
            align.mapq,
            align_score,
        )


@pd.api.extensions.register_dataframe_accessor("phased_bamdf")
class PhasedBamEntryDf(basePorecDf):
    """An extension to handle a dataframe derived from a BAM file"""

    DTYPE = {**BamEntryDf.DTYPE, **{"phase_block": GENOMIC_COORD_DTYPE, "haplotype": HAPLOTYPE_IDX_DTYPE}}

    @staticmethod
    def align_to_tuple(read_idx: int, align_idx: int, align: AlignedSegment) -> Tuple:
        t = BamEntryDf.align_to_tuple(read_idx, align_idx, align)
        try:
            haplotype = align.get_tag("HP")
            phase_block = align.get_tag("PS")
        except KeyError:
            haplotype = -1
            phase_block = 0
        res = (*t, *(phase_block, haplotype))
        return res


@pd.api.extensions.register_dataframe_accessor("porec_align")
class PoreCAlignDf(basePorecDf):
    """An extension to handle poreC-annotated alignments"""

    DTYPE = {
        "read_idx": READ_IDX_DTYPE,
        "align_idx": ALIGN_IDX_DTYPE,
        "mapping_type": pd.CategoricalDtype(["unmapped", "primary", "supplementary", "secondary"], ordered=True),
        "chrom": str,
        "start": GENOMIC_COORD_DTYPE,
        "end": GENOMIC_COORD_DTYPE,
        "strand": bool,
        "read_name": str,
        "read_length": READ_COORD_DTYPE,
        "read_start": READ_COORD_DTYPE,
        "read_end": READ_COORD_DTYPE,
        "mapping_quality": np.uint8,
        "score": np.uint32,
        # fields above come from BAM, below are calculated by pore-c tools
        "pass_filter": bool,
        "reason": pd.CategoricalDtype(
            ["pass", "unmapped", "singleton", "low_mq", "overlap_on_read", "not_on_shortest_path"], ordered=True
        ),
        "fragment_id": FRAG_IDX_DTYPE,
        "contained_fragments": np.uint32,
        "fragment_start": GENOMIC_COORD_DTYPE,
        "fragment_end": GENOMIC_COORD_DTYPE,
        "perc_of_alignment": PERCENTAGE_DTYPE,
        "perc_of_fragment": PERCENTAGE_DTYPE,
    }
    NANS = {
        "fragment_id": 0,
        "contained_fragments": 0,
        "fragment_start": 0,
        "fragment_end": 0,
        "perc_of_alignment": -1.0,
        "perc_of_fragment": -1.0,
    }


@pd.api.extensions.register_dataframe_accessor("phased_porec_align")
class PhasedPoreCAlignDf(basePorecDf):
    """An extension to handle poreC-annotated alignments"""

    DTYPE = {**PoreCAlignDf.DTYPE, **{"phase_block": GENOMIC_COORD_DTYPE, "haplotype": HAPLOTYPE_IDX_DTYPE}}
    NANS = {**PoreCAlignDf.NANS, **{"phase_block": -1, "haplotype": 0}}


@pd.api.extensions.register_dataframe_accessor("porec_read")
class PoreCReadDf(basePorecDf):
    DTYPE = {
        "read_idx": READ_IDX_DTYPE,
        "read_name": str,
        "read_length": READ_COORD_DTYPE,
        "num_aligns": np.uint16,
        "num_pass_aligns": np.uint16,
        "unique_fragments_assigned": np.uint16,
        "contained_fragments": np.uint16,
        "num_contacts": np.uint32,
        "num_cis_contacts": np.uint32,
        "perc_read_assigned": PERCENTAGE_DTYPE,
        "num_chroms_contacted": np.uint16,
    }
    NANS = {
        "num_aligns": 0,
        "num_pass_aligns": 0,
        "unique_fragments_assigned": 0,
        "contained_fragments": 0,
        "num_contacts": 0,
        "num_cis_contacts": 0,
        "perc_read_assigned": -1.0,
        "num_chroms_contacted": 0,
    }


@pd.api.extensions.register_dataframe_accessor("fragmentdf")
class FragmentDf(basePorecDf):
    """An extension to handle dataframes containing pairfile data"""

    DTYPE = {
        "fragment_id": FRAG_IDX_DTYPE,  # uid starting at 1
        "chrom": None,  # will be categorical but unknown until runtimei, use set_dtype
        "start": GENOMIC_COORD_DTYPE,
        "end": GENOMIC_COORD_DTYPE,
        "fragment_length": GENOMIC_COORD_DTYPE,
    }

    def assert_valid(self):
        if "chrom" in self._validation_errors and self._validation_errors["chrom"].startswith("Mismatched dtype"):
            self._validation_errors.pop("chrom")
        super(FragmentDf, self).assert_valid()


@pd.api.extensions.register_dataframe_accessor("pairdf")
class PairDf(object):
    """An extension to handle dataframes containing pairfile data"""

    DTYPE = {
        "readID": str,
        "chr1": str,
        "pos1": np.uint32,
        "chr2": str,
        "pos2": np.uint32,
        "strand1": pd.CategoricalDtype(["+", "-"]),
        "strand2": pd.CategoricalDtype(["+", "-"]),
        "pair_type": pd.CategoricalDtype(["DJ", "IJ"]),  # DJ=direct_junction, IJ=indirect_junction
        "frag1": np.uint32,
        "frag2": np.uint32,
        "align_idx1": np.uint32,
        "align_idx2": np.uint32,
        "distance_on_read": np.int32,
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert obj.dtype == PairDf.DTYPE

    def is_valid(self):
        # FIXFIX
        return True


@pd.api.extensions.register_dataframe_accessor("contactdf")
class ContactDf(object):
    """An extension to handle dataframes containing pairfile data"""

    DTYPE = {
        "read_name": str,
        "read_idx": READ_IDX_DTYPE,
        "read_length": READ_COORD_DTYPE,
        "read_order": np.uint16,
        "contact_is_direct": bool,
        "contact_is_cis": bool,
        "contact_genome_distance": GENOMIC_COORD_DTYPE,
        "contact_fragment_distance": GENOMIC_COORD_DTYPE,
        # align1
        "align1_align_idx": np.uint32,
        "align1_chrom": str,
        "align1_start": GENOMIC_COORD_DTYPE,
        "align1_end": GENOMIC_COORD_DTYPE,
        "align1_strand": bool,
        "align1_read_start": READ_COORD_DTYPE,
        "align1_read_end": READ_COORD_DTYPE,
        "align1_mapping_quality": np.uint8,
        "align1_score": np.uint32,
        "align1_fragment_id": FRAG_IDX_DTYPE,
        "align1_fragment_start": GENOMIC_COORD_DTYPE,
        "align1_fragment_end": GENOMIC_COORD_DTYPE,
        "align1_fragment_midpoint": GENOMIC_COORD_DTYPE,
        # align2
        "align2_align_idx": np.uint32,
        "align2_chrom": str,
        "align2_start": GENOMIC_COORD_DTYPE,
        "align2_end": GENOMIC_COORD_DTYPE,
        "align2_strand": bool,
        "align2_read_start": READ_COORD_DTYPE,
        "align2_read_end": READ_COORD_DTYPE,
        "align2_mapping_quality": np.uint8,
        "align2_score": np.uint32,
        "align2_fragment_id": FRAG_IDX_DTYPE,
        "align2_fragment_start": GENOMIC_COORD_DTYPE,
        "align2_fragment_end": GENOMIC_COORD_DTYPE,
        "align2_fragment_midpoint": GENOMIC_COORD_DTYPE,
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert obj.dtype == PairDf.DTYPE

    def is_valid(self):
        return True


@pd.api.extensions.register_dataframe_accessor("phased_contactdf")
class PhasedContactDf(object):
    """An extension to handle dataframes containing pairfile data"""

    DTYPE = {
        **ContactDf.DTYPE,
        **{
            "align1_phase_block": GENOMIC_COORD_DTYPE,
            "align1_haplotype": HAPLOTYPE_IDX_DTYPE,
            "align2_phase_block": GENOMIC_COORD_DTYPE,
            "align2_haplotype": HAPLOTYPE_IDX_DTYPE,
        },
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert obj.dtype == PairDf.DTYPE

    def is_valid(self):
        return True


@pd.api.extensions.register_dataframe_accessor("salsadf")
class SalsaDf(object):
    """An extension to handle dataframes containing salsa bed format data"""

    DTYPE = {
        "chr": str,
        "start": GENOMIC_COORD_DTYPE,
        "end": GENOMIC_COORD_DTYPE,
        "read_pair_id": str,
        "mapping_quality": np.uint8,
        "strand": pd.CategoricalDtype(["+", "-"]),
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert obj.dtype == PairDf.DTYPE

    def is_valid(self):
        return True


@pd.api.extensions.register_dataframe_accessor("hictxt")
class HicTxtDf(object):
    """An extension to handle dataframes containing .hic data"""

    DTYPE = {
        "readID": str,
        "strand1": pd.CategoricalDtype(["0", "16"]),
        "chr1": str,
        "pos1": GENOMIC_COORD_DTYPE,
        "frag1": np.uint32,
        "strand2": pd.CategoricalDtype(["0", "16"]),
        "chr2": str,
        "pos2": GENOMIC_COORD_DTYPE,
        "frag2": np.uint32,
        "mapping_quality1": np.uint8,
        "mapping_quality2": np.uint8,
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert obj.dtype == PairDf.DTYPE

    def is_valid(self):
        return True


@pd.api.extensions.register_dataframe_accessor("aligndf")
class AlignDf(object):
    """An extension to handle dataframes containing alignment intervals"""

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        for col in ["chrom", "start", "end"]:
            if col not in obj.columns:
                raise AttributeError("Must have columns 'chrom', 'start' and 'end'.")
        self.index_name = "index" if obj.index.name is None else obj.index.name
        assert obj.index.is_unique, "Must have a unique index"
        assert not isinstance(obj.index, pd.MultiIndex), "Can't be multindex"
        assert np.issubdtype(obj.index.dtype, np.integer), "Must have integer index: {}".format(obj.index.dtype)


@pd.api.extensions.register_dataframe_accessor("ginterval")
class GenomeIntervalDf(object):
    """An extension to handle dataframes containing genomic intervals"""

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        for col in ["chrom", "start", "end"]:
            if col not in obj.columns:
                raise AttributeError("Must have columns 'chrom', 'start' and 'end'.")
        self.index_name = "index" if obj.index.name is None else obj.index.name
        dupes = obj.index.duplicated(keep=False)
        if dupes.any():
            raise ValueError("Must have a unique index: {}".format(obj[dupes]))
        assert not isinstance(obj.index, pd.MultiIndex), "Can't be multindex"
        if not np.issubdtype(obj.index.dtype, np.integer):
            raise ValueError("Must have integer index: {}\n{}".format(obj.index.dtype, obj))

    @property
    def is_valid(self):
        return True

    def as_ncls_dict(self) -> Dict[Chrom, NCLS]:
        res = {}
        for chrom, chrom_df in self._obj.groupby("chrom"):
            res[chrom] = NCLS(
                chrom_df.start.values.astype(np.int64),
                chrom_df.end.values.astype(np.int64),
                chrom_df.index.values.astype(np.int64),
            )
        return res

    def assign(self, other: "GenomeIntervalDf"):
        """Assign intervals in this dataframe to the first overlapping interval in other"""
        tgt = other.ginterval.as_ncls_dict()
        overlaps = []
        for chrom, chrom_df in self._obj.groupby("chrom"):
            if chrom not in tgt:
                continue
            self_indices, target_indices = tgt[chrom].first_overlap_both(
                chrom_df.start.values.astype(np.int64),
                chrom_df.end.values.astype(np.int64),
                chrom_df.index.values.astype(np.int64),
            )
            if len(self_indices) == 0:
                continue
            overlaps.append(pd.DataFrame({"self": self_indices, "other": target_indices}))
        overlaps = pd.concat(overlaps, ignore_index=True).set_index("self")
        if overlaps.index.duplicated().any():
            raise ValueError(overlaps[overlaps.index.duplicated(keep="both")])
        return overlaps.reindex(index=self._obj.index)

    def overlap(self, other: "GenomeIntervalDf", calculate_lengths: bool = True):
        """Find all overlapping intervals between this dataframe and 'other'"""
        self_cols = ["chrom", "start", "end", self.index_name]
        other_rename = {"start": "other_start", "end": "other_end"}
        if self.index_name == other.ginterval.index_name:
            other_rename[other.ginterval.index_name] = "other_" + other.ginterval.index_name
        else:
            other_rename[other.ginterval.index_name] = other.ginterval.index_name
        tgt = other.ginterval.as_ncls_dict()
        overlaps = []
        for chrom, chrom_df in self._obj.groupby("chrom"):
            if chrom not in tgt:
                continue
            self_indices, target_indices = tgt[chrom].all_overlaps_both(
                chrom_df.start.values.astype(np.int64),
                chrom_df.end.values.astype(np.int64),
                chrom_df.index.values.astype(np.int64),
            )
            overlap_df = (
                chrom_df.reindex(self_indices).reset_index()
                # .rename(columns={'index': 'query_idx'})
                .loc[:, self_cols]
            ).join(
                other.reindex(target_indices).reset_index().rename(columns=other_rename).loc[:, other_rename.values()]
            )
            if calculate_lengths:
                overlap_df = (
                    overlap_df.assign(
                        overlap_start=lambda x: np.where(x.other_start > x.start, x.other_start, x.start).astype(int),
                        overlap_end=lambda x: np.where(x.other_end < x.end, x.other_end, x.end).astype(int),
                    )
                    .eval("overlap_length = overlap_end - overlap_start")
                    .eval("perc_of_self = (100.0 * overlap_length) / (end - start)")
                    .eval("perc_of_other = (100.0 * overlap_length) / (other_end - other_start)")
                )
            overlaps.append(overlap_df)
        if overlaps:
            res = pd.concat(overlaps, ignore_index=True).astype({self.index_name: np.uint64})
        else:
            res = None
        return res

    @classmethod
    def fixed_width_bins(cls, chrom_lengths, bin_width):
        dfs = []
        for chrom_id, chrom_length in chrom_lengths.items():
            bin_edges = list(range(0, chrom_length, bin_width))
            if bin_edges[-1] != chrom_length:
                bin_edges.append(chrom_length)
            _df = (
                pd.DataFrame({"start": bin_edges[:-1], "end": bin_edges[1:]})
                .astype(GENOMIC_COORD_DTYPE)
                .assign(chrom=chrom_id)
            )
            dfs.append(_df)
        df = pd.concat(dfs, ignore_index=True).reset_index().rename(columns={"index": "bin_id"})
        return df[["chrom", "start", "end", "bin_id"]]
