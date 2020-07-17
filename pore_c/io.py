import subprocess as sp
from logging import getLogger
from pathlib import Path

import pyarrow as pa
from pyarrow import parquet as pq


logger = getLogger(__name__)


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


class BatchedFastqWriter(object):
    def __init__(self, output_path: Path):
        self._output_path = output_path
        if self._output_path.suffix == ".gz":
            self._raw_output_path = self._output_path.with_suffix("")
            self._compress = True
        else:
            self._raw_output_path = self._output_path
            self._compress = False

        self.output_paths = []
        self._counter = 0
        self._batch_counter = 0

    def __call__(self, sequences):
        if len(sequences) == 0:
            return
        self._batch_counter += 1
        output_path = str(self._raw_output_path).format(self._batch_counter)
        with open(output_path, "w") as fh:
            fh.write("%s\n" % "\n".join(sequences))
        if self._compress:
            logger.info("Compressing {} using bgzip".format(output_path))
            sp.check_call(["bgzip", str(output_path)])
        self.output_paths.append(output_path)
        self._counter += len(sequences)

    def close(self):
        logger.debug(
            "Wrote {} sequences to {} files: {}".format(self._counter, self._batch_counter, " ".join(self.output_paths))
        )


class TableWriter(object):
    def __init__(self, path, version="2.0"):
        self.path = path
        self.writer = None
        self.schema = None
        self.counter = 0
        self.row_counter = 0
        self.version = version

    def write(self, df):
        table = pa.Table.from_pandas(df, preserve_index=False, schema=self.schema)
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.path, schema=table.schema, version=self.version)
            self.schema = table.schema
        try:
            self.writer.write_table(table)
        except Exception as exc:
            raise IOError("Error writing batch {} to {}:\n{}\n{}".format(self.counter, self.path, df.head(), exc))
        self.row_counter += len(df)
        self.counter += 1

    def __call__(self, *args, **kwds):
        return self.write(*args, **kwds)

    def close(self):
        self.writer.close()
