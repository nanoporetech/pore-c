import numpy as np
import pandas as pd
import pysam

from pore_c.analyses.reads import mean_qscore


def test_mean_qscore(read_fastq_table, read_fastq_file):
    """Check the mean qscore against that produced by seqkit"""
    expected = pd.read_csv(
        read_fastq_table, sep="\t", usecols=["read_id", "mean_qscore"], index_col="read_id", squeeze=True
    )
    for rec in pysam.FastqFile(read_fastq_file):
        qscore = mean_qscore(rec.get_quality_array())
        assert np.around(qscore, 2) == expected[rec.name]
