import sys
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd
from streamz import Stream

from pore_c.datasources import Fastq
from pore_c.io import FastqWriter
from pore_c.utils import DataFrameProgress, mean_qscore


logger = getLogger(__name__)


def read_length_stats(lengths, percentiles=[25, 50, 75]):
    if len(lengths) == 0:
        return {}
    sorted_lengths = np.sort(lengths.values)
    cumsum = sorted_lengths.cumsum()
    total_bases = cumsum[-1]
    n50_idx = np.argmax(cumsum > (total_bases * 0.5))
    n50 = sorted_lengths[n50_idx]
    qdict = dict(zip(["Q%d" % p for p in percentiles], map(float, np.percentile(sorted_lengths, percentiles))))
    return {
        **{
            "num_sequences": len(lengths),
            "total_bases": int(total_bases),
            "mean": float(total_bases / len(lengths)),
            "min": int(sorted_lengths[0]),
            "max": int(sorted_lengths[-1]),
            "N50": int(n50),
            **qdict,
        }
    }


def prepare_fastq(
    input_fastq: Path,
    pass_fastq: Path = None,
    fail_fastq: Path = None,
    read_metadata: Path = None,
    summary: Path = None,
    min_read_length: int = 50,
    max_read_length: int = 5000000,
    min_qscore: int = 0,
    max_qscore: int = 266,
    chunksize: int = 10000,
):

    fastq_stream = Stream()

    filtered_stream = fastq_stream.map(filter_records, min_read_length, max_read_length, min_qscore, max_qscore)

    pass_writer = FastqWriter(pass_fastq)
    fail_writer = FastqWriter(fail_fastq)

    read_prog = ReadFilterProgress()

    df_sink = (
        filtered_stream.pluck("metadata").accumulate(read_prog, returns_state=True, start=read_prog).sink_to_list()
    )
    pass_sink = filtered_stream.pluck("pass").sink(pass_writer)  # noqa: F841
    fail_sink = filtered_stream.pluck("fail").sink(fail_writer)  # noqa: F841

    # split reads into chunks for processing
    for chunk_idx, records in enumerate(Fastq(input_fastq).read_chunked(chunksize)):
        fastq_stream.emit(records)
    metadata_df = pd.concat(df_sink, ignore_index=True)
    metadata_df.to_parquet(read_metadata, index=False)
    pass_rate = metadata_df["pass_filter"].mean()
    if pass_rate == 0:
        raise ValueError("No reads passed filter")
    summary_stats = {
        "all": read_length_stats(metadata_df["read_length"]),
        "pass": read_length_stats(metadata_df.query("pass_filter == True")["read_length"]),
        "fail": read_length_stats(metadata_df.query("pass_filter == False")["read_length"]),
    }
    pass_writer.close()
    fail_writer.close()
    read_prog.close()
    sys.stderr.write("\n")
    logger.debug("Finished processing reads")
    assert (pass_writer._counter + fail_writer._counter) == summary_stats["all"]["num_sequences"]
    df = pd.DataFrame([v for v in summary_stats.values() if v], index=[k for k, v in summary_stats.items() if v])
    df.index.name = "read_subset"
    logger.info("Finished processing {}:\n{}\n".format(input_fastq, str(df)))
    df.to_csv(summary)
    return summary_stats


def filter_records(list_of_records, min_read_length, max_read_length, min_qscore, max_qscore):

    df = (
        pd.DataFrame(
            [(_.name, len(_.sequence), mean_qscore(_.get_quality_array())) for _ in list_of_records],
            columns=["read_id", "read_length", "qscore"],
        )
        .astype({"read_length": pd.np.uint32, "qscore": pd.np.float32})
        .eval(
            "pass_filter = (@min_read_length <= read_length < @max_read_length) & (@min_qscore <= qscore < @max_qscore)"
        )
    )
    seq_strings = {True: [], False: []}
    for seq, is_pass in zip(map(str, list_of_records), df["pass_filter"].values):
        seq_strings[is_pass].append(seq)

    return {"metadata": df, "pass": seq_strings[True], "fail": seq_strings[False]}


class ReadFilterProgress(DataFrameProgress):
    def __init__(self, **kwds):
        kwds["desc"] = "Reads processed"
        kwds["unit"] = " reads"
        super(ReadFilterProgress, self).__init__(**kwds)

    def update_data(self, read_df):
        count_cols = ["read_length", "pass_filter"]
        counts = read_df.loc[:, count_cols].sum(axis=0)
        counts["reads"] = len(read_df)
        counts = counts.astype(int).to_frame().T

        if self._data is None:
            self._data = counts
        else:
            self._data += counts

    @staticmethod
    def summarize(df):
        summary = (
            df.eval("Gb = read_length * 1e-9").eval("perc_pass = 100 * (pass_filter / reads)")
            # .loc[:, ["reads", "Gb", "num_contacts", "contacts_per_Gb", "percent_cis"]]
        )
        return summary

    def update_progress_bar(self, num_aligns):
        self._bar.update(num_aligns)
        self._bar.set_postfix(self.summarize(self._data).iloc[0, :].to_dict())
