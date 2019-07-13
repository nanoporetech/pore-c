from collections import defaultdict
from dataclasses import dataclass
from typing import (Dict, Iterator, List, NamedTuple, NewType, Pattern, Set,
                    Tuple, Union)

import numpy as np
import pandas as pd
from ncls import NCLS
from pybedtools import BedTool

Chrom = NewType("Chrom", str)



#int8 	Byte (-128 to 127)
#int16 	Integer (-32768 to 32767)
#int32 	Integer (-2147483648 to 2147483647)
#int64 	Integer (-9223372036854775808 to 9223372036854775807)
#uint8 	Unsigned integer (0 to 255)
#uint16 	Unsigned integer (0 to 65535)
#uint32 	Unsigned integer (0 to 4294967295)
#uint64 	Unsigned integer (0 to 18446744073709551615)


GENOMIC_COORD_DTYPE = np.uint32  # should be fine as long as individual chromosomes are less than 4Gb
READ_COORD_DTYPE = np.uint32
STRAND_DTYPE = pd.CategoricalDtype(['+','-'], ordered=True)
FRAG_IDX_DTYPE = np.uint32
READ_IDX_DTYPE = np.uint32
ALIGN_IDX_DTYPE = np.uint32


class basePorecDf(object):
    DTYPE = {}
    NANS = {}
    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        self._validation_errors = {}
        self._additional_columns = []
        self._check(pandas_obj)

    def _check(self, obj):
        errors = {}
        for key, value in self.DTYPE.items():
            _dtype = obj.dtypes.get(key, None)
            if _dtype is None:
                errors[key] = 'Missing column'
            elif value is None:
                errors[key] = 'Unset datatype'
            elif _dtype != value:
                errors[key] = "Mismatched dtype ({}/{}) (expected/found)".format(_dtype, value)
        self._validation_errors = errors
        for column, dtype in obj.dtypes.items():
            if column not in self.DTYPE:
                self._additional_columns.append(column)

    def assert_valid(self):
        if self._validation_errors:
            raise ValueError("Failed validation:\n{}\n".format('\n'.join(["{}: {}".format(*_) for _ in self._validation_errors.items()])))

    @classmethod
    def set_dtype(cls, key, value):
        """If we don't know a datatype until execution time we can set it here"""
        cls.DTYPE[key] = value

    def cast(self, subset=False, fillna=False) -> pd.DataFrame:
        cols = list(self.DTYPE.keys()) if subset else self._obj.columns
        fillna_dict = self.NANS if fillna else {}
        return self._obj.loc[:, cols].fillna(fillna_dict).astype(self.DTYPE)


@pd.api.extensions.register_dataframe_accessor("fragmentdf")
class FragmentDf(basePorecDf):
    """An extension to handle dataframes containing pairfile data"""
    DTYPE = {
        'fragment_id': FRAG_IDX_DTYPE, # uid starting at 1
        'chrom': None, # will be categorical but unknown until runtimei, use set_dtype
        'start': GENOMIC_COORD_DTYPE,
        'end': GENOMIC_COORD_DTYPE,
        'fragment_length': GENOMIC_COORD_DTYPE,
    }
    def assert_valid(self):
        if 'chrom' in self._validation_errors and self._validation_errors['chrom'].startswith('Mismatched dtype'):
            self._validation_errors.pop('chrom')
        super(FragmentDf, self).assert_valid()






@pd.api.extensions.register_dataframe_accessor("pairdf")
class PairDf(object):
    """An extension to handle dataframes containing pairfile data"""
    DTYPE = {
        'readID': str,
        'chr1': str,
        'pos1': np.uint32,
        'chr2': str,
        'pos2': np.uint32,
        'strand1': pd.CategoricalDtype(['+','-']),
        'strand2': pd.CategoricalDtype(['+','-']),
        'pair_type': pd.CategoricalDtype(['DJ', 'IJ']), # DJ=direct_junction, IJ=indirect_junction
        'frag1': np.uint32,
        'frag2': np.uint32,
        'align_idx1': np.uint32,
        'align_idx2': np.uint32,
        'distance_on_read': np.int32,
    }

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    def _validate(self, obj):
        assert(obj.dtype == PairDf.DTYPE)

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
        assert obj.index.is_unique, "Must have a unique index"
        assert not isinstance(obj.index, pd.MultiIndex), "Can't be multindex"
        assert np.issubdtype(obj.index.dtype, np.integer), "Must have integer index: {}".format(obj.index.dtype)

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

    def overlap(self, other: "GenomeIntervalDf", calculate_lengths: bool = True):
        """Find all overlapping intervals between this dataframe and 'other'"""
        self_cols = ["chrom", "start", "end", self.index_name]
        other_cols = ["start", "end", other.ginterval.index_name]
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


class SeqDigest(NamedTuple):
    """Holds the results of a digestion"""

    seq_name: str
    positions: List[int]
    seq_length: int


@dataclass
class BedToolsOverlap(object):
    """Datastructure to hold the results of a bedtools intersect command"""

    query_chrom: str
    query_start: int
    query_end: int
    query_id: str
    frag_chrom: str
    frag_start: int
    frag_end: int
    frag_id: int
    overlap: int

    def __post_init__(self):
        self.query_start = int(self.query_start)
        self.query_end = int(self.query_end)
        self.frag_start = int(self.frag_start)
        self.frag_end = int(self.frag_end)
        self.frag_id = int(self.frag_id)
        self.overlap = int(self.overlap)


class FragmentMap(object):
    """Represents the fragments created by a restriction digestion"""

    def __init__(self, bedtool: BedTool, chrom_lengths: Dict[str, int] = None, Terminal_Fragments: Set = None):
        self.bt = bedtool  # bedtool.saveas().tabix(is_sorted=True)
        self.chrom_lengths = chrom_lengths

        # this will store a list of fragment_ids which, while adjacent in number space,
        # are not physically adjacent, e.g. q-telomere of chr1 and p-telomere of chr2.
        self.terminal_fragments = Terminal_Fragments

    @staticmethod
    def endpoints_to_intervals(chrom, positions, id_offset, chrom_length=None) -> List[Tuple[str, int, int, int]]:
        if len(positions) == 0:  # if there's no cut sites on a sequence then we return the whole sequence
            assert chrom_length is not None
            positions = [0, chrom_length]
        else:
            if not positions[0] == 0:
                positions = [0] + positions
            if chrom_length and positions[-1] != chrom_length:
                positions = positions + [chrom_length]
        return [
            (chrom, start, end, str(x + id_offset)) for x, (start, end) in enumerate(zip(positions[:-1], positions[1:]))
        ]

    def save_to_bed(self, path):
        if not path.endswith(".bed.gz"):
            raise ValueError("Must end with .bed.gz: {}".format(path))
        self.bt.saveas(path.rsplit(".gz")[0]).tabix(in_place=True, is_sorted=True, force=True)

    def save_to_HiCRef(self, path):
        chrom_to_endpoints = defaultdict(list)
        for interval in self.bt:
            chrom_to_endpoints[interval.chrom].append(interval.stop)
        with open(path, "w") as fh:
            for chrom, endpoints in chrom_to_endpoints.items():
                fh.write("{} {}\n".format(chrom, " ".join(map(str, endpoints))))

    @classmethod
    def from_dict(cls, d: Dict[str, List[int]]):
        """Create a fragment map from a dictionary mapping chromosomes to fragment endpoints (useful for testing)"""
        intervals = []
        chrom_lengths = {}
        id_offset = 0
        for chrom, endpoints in d.items():
            intervals.extend(cls.endpoints_to_intervals(chrom, endpoints, id_offset))
            chrom_lengths[chrom] = endpoints[-1]
            id_offset += len(endpoints)
        bt = BedTool(intervals)
        return cls(bt, chrom_lengths)

    @classmethod
    def from_bed_file(cls, path):
        if not path.endswith(".bed.gz"):
            raise ValueError("Must end with .bed.gz: {}".format(path))
        return cls(BedTool(path))

    @classmethod
    def from_HiCRef(cls, fname):
        intervals = []
        chrom_lengths = {}
        id_offset = 0
        with open(fname) as fh:
            for line in fh:
                fields = line.strip().split()
                chrom = fields[0]
                endpoints = list(map(int, fields[1:]))
                intervals.extend(cls.endpoints_to_intervals(chrom, endpoints, id_offset))
                chrom_lengths[chrom] = endpoints[-1]
                id_offset += len(endpoints)
        bt = BedTool(intervals)
        return cls(bt, chrom_lengths, terminal_fragments)

    @classmethod
    def from_digest_iter(cls, i: Iterator[SeqDigest]):
        intervals = []
        chrom_lengths = {}
        terminal_fragments = set()
        id_offset = 0
        for digest in i:
            chrom = digest.seq_name
            endpoints = digest.positions
            new_intervals = cls.endpoints_to_intervals(chrom, endpoints, id_offset, chrom_length=digest.seq_length)
            intervals.extend(new_intervals)
            chrom_lengths[chrom] = digest.seq_length
            id_offset += len(endpoints)
            if id_offset != 0:
                terminal_fragments.add((id_offset, id_offset + 1))

        bt = BedTool(intervals)
        return cls(bt, chrom_lengths, terminal_fragments)

    def _query_to_bedtool(self, query):
        def _interval_from_tuple(t, id_offset=0):
            if len(t) == 3:
                return (t[0], t[1], t[2], str(id_offset))
            else:
                assert len(t) == 4
                return t

        if isinstance(query, tuple):
            intervals = [_interval_from_tuple(query)]
        else:
            raise NotImplementedError
        return BedTool(intervals)

    def iter_overlaps(self, query, min_overlap=0):
        query_bt = self._query_to_bedtool(query)
        for overlap in (BedToolsOverlap(*i.fields) for i in query_bt.intersect(self.bt, sorted=True, wo=True)):
            if overlap.overlap <= min_overlap:
                continue
            yield overlap
