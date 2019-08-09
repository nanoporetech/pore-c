from typing import Dict, NewType

import numpy as np
import pandas as pd
from ncls import NCLS

Chrom = NewType("Chrom", str)


# int8 	Byte (-128 to 127)
# int16 	Integer (-32768 to 32767)
# int32 	Integer (-2147483648 to 2147483647)
# int64 	Integer (-9223372036854775808 to 9223372036854775807)
# uint8 	Unsigned integer (0 to 255)
# uint16 	Unsigned integer (0 to 65535)
# uint32 	Unsigned integer (0 to 4294967295)
# uint64 	Unsigned integer (0 to 18446744073709551615)
GENOMIC_COORD_DTYPE = np.uint32  # should be fine as long as individual chromosomes are less than 4Gb
READ_COORD_DTYPE = np.uint32
STRAND_DTYPE = pd.CategoricalDtype(["+", "-"], ordered=True)
FRAG_IDX_DTYPE = np.uint32
READ_IDX_DTYPE = np.uint32
ALIGN_IDX_DTYPE = np.uint32
PERCENTAGE_DTYPE = np.float32


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
