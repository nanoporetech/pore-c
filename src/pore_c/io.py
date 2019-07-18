import subprocess as sp
from logging import getLogger

import pyarrow as pa
from pyarrow import parquet as pq

logger = getLogger(__name__)

class PairFileWriter(object):
    def __init__(self, output_path, chrom_sizes, genome_assembly, columns=None):
        self._output_path = output_path
        if self._output_path.suffix == ".gz":
            self._raw_output_path = self._output_path.with_suffix("")
            self._sort_and_compress = True
        else:
            self._raw_output_path = self._output_path
            self._sort_and_compress = False
        self._chrom_sizes = chrom_sizes
        self._genome_assembly = genome_assembly
        self._columns = columns
        self._fh = None
        self._string_template = "%s\n" % ("\t".join(["{}" for c in self._columns]))

    def _write_header(self):
        lines = [
            "## pairs format 1.0",
            "#shape: upper triangle",
            "#genome_assembly {}".format(self._genome_assembly),
        ]
        lines.extend(
            ["#chromsize: {} {}".format(*_) for _ in self._chrom_sizes.items()]
        )
        lines.append("#columns: {}".format(" ".join(self._columns)))
        self._fh.write("{}\n".format("\n".join(lines)))

    def __call__(self, pair_df):
        assert len(pair_df.columns) == len(self._columns)
        if self._fh is None:
            self._fh = open(self._raw_output_path, "w")
            self._write_header()
        pair_df.to_csv(self._fh, header=None, sep="\t", index=False)

    def close(self):
        self._fh.close()
        # TODO: this could all be cleaned up a bit
        if self._sort_and_compress:
            logger.info("Sorting and compressing: {}".format(self._output_path))
            comd = "pairtools sort {} | bgzip > {}".format(
                self._raw_output_path, self._output_path
            )
            logger.debug('Running command: {}'.format(comd))
            sp.check_call(comd, shell=True)
            sp.check_call(["pairix", str(self._output_path)])
            logger.debug('Removing temp file: {}'.format(self._raw_output_path))
            sp.check_call(["rm", str(self._raw_output_path)])


class FastqWriter(object):
    def __init__(self, output_path):
        self._output_path = output_path
        if self._output_path.suffix == ".gz":
            self._raw_output_path = self._output_path.with_suffix("")
            self._compress = True
        else:
            self._raw_output_path = self._output_path
            self._compress = False
        self._fh = None
        self._counter = 0

    def __call__(self, sequences):
        if len(sequences) == 0:
            return
        if self._fh is None:
            self._fh = open(self._raw_output_path, "w")
        self._fh.write("%s\n" % "\n".join(sequences))
        self._counter += len(sequences)

    def close(self):
        if self._fh is None:
            assert self._counter == 0
            return
        self._fh.close()
        logger.debug("Wrote {} sequences to {}".format(self._counter, self._raw_output_path))
        # TODO: this could all be cleaned up a bit
        if self._compress:
            logger.info("Compressing {} using bgzip".format(self._raw_output_path))
            sp.check_call(["bgzip", str(self._raw_output_path)])


class TableWriter(object):
    def __init__(self, path):
        self.path = path
        self._writer = None
        self._schema = None
        self._counter = 0

    def write(self, df):
        table = pa.Table.from_pandas(df, preserve_index=False, schema=self._schema)
        if self._writer is None:
            self._writer = pq.ParquetWriter(self.path, schema=table.schema)
            self._schema = table.schema
        try:
            self._writer.write_table(table)
        except Exception as exc:
            raise IOError(
                "Error writing batch {} to {}:\n{}\n{}".format(
                    self._counter, self.path, df.head(), exc
                )
            )
        self._counter += 1

    def __call__(self, *args, **kwds):
        return self.write(*args, **kwds)

    def close(self):
        self._writer.close()
