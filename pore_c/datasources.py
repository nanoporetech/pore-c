import dask
from intake.source.base import DataSource, Schema
from pysam import FastaFile, FastxFile


class Fastq(DataSource):
    name = "fastq"
    version = "0.1.0"
    container = "python"
    partition_access = False
    description = "A fastq file"

    def __init__(self, urlpath, metadata=None):
        self._urlpath = urlpath
        super().__init__(metadata=metadata)

    def _open_dataset(self):
        return FastxFile(self._urlpath)

    def _get_schema(self):
        return Schema(datashape=None, dtype=None, shape=None, npartitions=None, extra_metadata={})

    def read_chunked(self, chunksize=10000):
        self._load_metadata()
        from toolz import partition_all

        yield from partition_all(chunksize, self._open_dataset())

    def _close(self):
        if self._dataset is not None:
            self._dataset.close()


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
        super().__init__(metadata=metadata)

    def _open_dataset(self):
        self._dataset = FastaFile(self._urlpath)

    def _get_schema(self):
        if self._dataset is None:
            self._open_dataset()
        self._chroms = list(self._dataset.references)
        chrom_lengths = [{"chrom": t[0], "length": t[1]} for t in zip(self._dataset.references, self._dataset.lengths)]
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
        return db.from_delayed([dask.delayed(self._get_partition(i)) for i in range(self.npartitions)])

    def _close(self):
        # close any files, sockets, etc
        if self._dataset is not None:
            self._dataset.close()
