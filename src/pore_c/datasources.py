from typing import Iterator, List, Tuple

import dask
import numpy as np
import pandas as pd
from intake.source.base import DataSource, Schema
from pysam import AlignedSegment, AlignmentFile, FastaFile, TabixFile, asTuple

from pore_c.model import BamEntryDf


class IndexedFasta(DataSource):
    name = "indexed_bedfile"
    version = "0.1.0"
    container = "python"
    partition_access = True
    description = "A bgzipped and indexed fasta file"

    def __init__(self, urlpath, metadata=None):
        self._urlpath = urlpath
        self._dataset = None
        self._dtype = None
        self._chroms = None
        super(IndexedFasta, self).__init__(metadata=metadata)

    def _open_dataset(self):
        self._dataset = FastaFile(self._urlpath)

    def _get_schema(self):
        if self._dataset is None:
            self._open_dataset()
        self._chroms = list(self._dataset.references)
        chrom_lengths = [
            {"chrom": t[0], "length": t[1]}
            for t in zip(self._dataset.references, self._dataset.lengths)
        ]
        return Schema(
            datashape=None,
            dtype=None,
            shape=None,
            npartitions=len(self._chroms),
            extra_metadata={"chroms": chrom_lengths},
        )

    def _get_partition(self, i):
        chrom = self._chroms[i]
        return [{"seqid": chrom, "seq": self._dataset.fetch(chrom)}]

    def read_chunked(self):
        self._load_metadata()
        for i in range(self.npartitions):
            yield self._get_partition(i)

    def to_dask(self):
        from dask import bag as db

        self._load_metadata()
        return db.from_delayed(
            [dask.delayed(self._get_partition(i)) for i in range(self.npartitions)]
        )

    def _close(self):
        # close any files, sockets, etc
        if self._dataset is not None:
            self._dataset.close()


class IndexedPairFile(DataSource):
    name = "indexed_pairfile"
    version = "0.1.0"
    container = "dataframe"
    partition_access = False
    description = "A bgzipped and indexed pairfile"

    def __init__(self, urlpath, include_unmapped=True, metadata=None):
        self._urlpath = urlpath
        self._include_unmapped = include_unmapped
        self._dataset = None
        self._dtype = None
        self._chroms = None
        super(IndexedPairFile, self).__init__(metadata=metadata)

    def _open_dataset(self):
        raise NotImplementedError
        self._dataset = TabixFile(self._urlpath)

    def _get_schema(self):
        raise NotImplementedError
        if self._dataset is None:
            self._open_dataset()
        self._chroms = list(self._dataset.contigs)

        rec = next(self._dataset.fetch(self._chroms[0], parser=asTuple()))
        num_fields = len(rec)

        chrom_coord_dtype = np.int64
        dtypes = {
            "chrom": pd.CategorialDtype(self._chroms + ["NULL"], ordered=True),
            "start": chrom_coord_dtype,
            "end": chrom_coord_dtype,
            "name": str,
            "score": np.float32,
            "strand": bool,
        }
        self._dtype = {key: dtypes[key] for key in list(dtypes.keys())[:num_fields]}
        return Schema(
            datashape=None,
            dtype=self._dtype,
            shape=(None, len(self._dtype)),
            npartitions=len(self._chroms),
            extra_metadata={},
        )

    def _get_partition(self, i):
        raise NotImplementedError
        chrom = self._chroms[i]
        columns = list(self._dtype.keys())
        return pd.DataFrame(
            list(self._dataset.fetch(chrom, parser=asTuple())), columns=columns
        ).astype(self._dtype)

    def read(self):
        raise NotImplementedError
        self._load_metadata()
        return pd.concat(
            [self.read_partition(i) for i in range(self.npartitions)], ignore_index=True
        )

    def _close(self):
        # close any files, sockets, etc
        if self._dataset is not None:
            self._dataset.close()


class IndexedBedFile(DataSource):
    name = "indexed_bedfile"
    version = "0.1.0"
    container = "dataframe"
    partition_access = False
    description = "A bgzipped and indexed bedfile"

    def __init__(self, urlpath, include_unmapped=True, metadata=None):
        self._urlpath = urlpath
        self._include_unmapped = include_unmapped
        self._dataset = None
        self._dtype = None
        self._chroms = None
        super(IndexedBedFile, self).__init__(metadata=metadata)

    def _open_dataset(self):
        self._dataset = TabixFile(self._urlpath)

    def _get_schema(self):
        if self._dataset is None:
            self._open_dataset()
        self._chroms = list(self._dataset.contigs)

        rec = next(self._dataset.fetch(self._chroms[0], parser=asTuple()))
        num_fields = len(rec)

        chrom_coord_dtype = np.int64
        dtypes = {
            "chrom": pd.CategorialDtype(self._chroms + ["NULL"], ordered=True),
            "start": chrom_coord_dtype,
            "end": chrom_coord_dtype,
            "name": str,
            "score": np.float32,
            "strand": bool,
        }
        self._dtype = {key: dtypes[key] for key in list(dtypes.keys())[:num_fields]}
        return Schema(
            datashape=None,
            dtype=self._dtype,
            shape=(None, len(self._dtype)),
            npartitions=len(self._chroms),
            extra_metadata={},
        )

    def _get_partition(self, i):
        chrom = self._chroms[i]
        columns = list(self._dtype.keys())
        return pd.DataFrame(
            list(self._dataset.fetch(chrom, parser=asTuple())), columns=columns
        ).astype(self._dtype)

    def read(self):
        self._load_metadata()
        return pd.concat(
            [self.read_partition(i) for i in range(self.npartitions)], ignore_index=True
        )

    def _close(self):
        # close any files, sockets, etc
        if self._dataset is not None:
            self._dataset.close()


class NameSortedBamSource(DataSource):
    name = "name_sorted_bam"
    version = "0.1.0"
    container = "dataframe"
    partition_access = False
    description = "Readname-sorted BAM of poreC alignments"

    def __init__(self, urlpath, include_unmapped=True, metadata=None):
        self._urlpath = urlpath
        self._include_unmapped = include_unmapped
        self._af = None
        self._dtype = None
        super(NameSortedBamSource, self).__init__(metadata=metadata)

    def _open_dataset(self):
        # TODO check that bam file has namesorted in header
        self._af = AlignmentFile(self._urlpath)

    def _get_schema(self):
        if self._af is None:
            self._open_dataset()
        chrom_names = list(self._af.references)
        assert "NULL" not in chrom_names
        dtype = BamEntryDf.DTYPE.copy()
        dtype["chrom"] = pd.CategoricalDtype(chrom_names + ["NULL"], ordered=True)
        self._dtype = dtype
        return Schema(
            datashape=None,
            dtype=dtype,
            shape=(None, len(dtype)),
            npartitions=None,
            extra_metadata={},
        )

    def get_chrom_dtype(self):
        return self._schema.dtype["chrom"]

    @staticmethod
    def _group_by_read(
        align_iter: Iterator[AlignedSegment]
    ) -> Iterator[List[Tuple[int, int, AlignedSegment]]]:
        """Iterate over alignments in name-sorted bam file grouping by read"""
        current_read_name = None
        read_idx = 0
        aligns = []
        for align_idx, align in enumerate(align_iter):
            if current_read_name is None:
                current_read_name = align.query_name
                aligns.append((read_idx, align_idx, align))
            elif current_read_name == align.query_name:
                aligns.append((read_idx, align_idx, align))
            else:
                yield aligns
                read_idx += 1
                current_read_name = align.query_name
                aligns = [(read_idx, align_idx, align)]
        yield aligns

    def _align_to_tuple(self, align_data):
        read_idx, align_idx, align = align_data
        if align.is_unmapped:
            align_cat = "unmapped"
            chrom, start, end, align_score = "NULL", 0, 0, 0
            read_length = align.query_length
        else:
            chrom, start, end = (
                align.reference_name,
                align.reference_start,
                align.reference_end,
            )
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
            align.query_name,
            read_length,
            align.query_alignment_start,
            align.query_alignment_end,
            align.mapq,
            align_score,
        )

    def read_chunked(self, chunksize=10000, yield_aligns=False, max_chunks=None):
        """Read the bam into a dataframe in chunks, yield those chunks
        """
        self._load_metadata()
        from toolz import partition_all

        BamEntryDf.set_dtype("chrom", self.get_chrom_dtype())
        align_iter = self._af.fetch(until_eof=self._include_unmapped)
        columns = list(self._schema.dtype.keys())
        for chunk_idx, chunk in enumerate(
            partition_all(chunksize, self._group_by_read(align_iter))
        ):
            aligns = [a for read_aligns in chunk for a in read_aligns]
            df = (
                pd.DataFrame([self._align_to_tuple(a) for a in aligns], columns=columns)
                .astype(BamEntryDf.DTYPE)
                .bamdf.cast(fillna=True, subset=True)
            )
            if yield_aligns:
                yield (aligns, df)
            else:
                yield (df)
            if max_chunks and chunk_idx == max_chunks - 1:
                break

    def _close(self):
        if self._af is not None:
            self._af.close()
