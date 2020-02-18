import logging
from collections import defaultdict

import dask.dataframe as dd
import numpy as np
import pandas as pd
import pyranges as pr
from cooler import Cooler, create_cooler
from cooler.util import binnify
from pysam import FastaFile

from pore_c.analyses.reference import revcomp
from pore_c.config import PQ_ENGINE, PQ_VERSION


logger = logging.getLogger(__name__)


def export_to_cooler(
    contact_table, output_prefix, cooler_resolution, fragment_table, chromsizes, query, query_columns=None
):
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []

    columns.extend(["align1_fragment_id", "align2_fragment_id"])
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns)
    if query:
        contact_df = contact_df.query(query)
    cooler_path = output_prefix + ".cool"
    chrom_dict = pd.read_csv(
        chromsizes, sep="\t", header=None, names=["chrom", "size"], index_col=["chrom"], squeeze=True
    )
    bins_df = binnify(chrom_dict, cooler_resolution)
    bins_df.index.name = "bin_id"
    bins = pr.PyRanges(bins_df.reset_index().rename(columns={"start": "Start", "end": "End", "chrom": "Chromosome"}))
    fragment_df = dd.read_parquet(fragment_table, engine=PQ_ENGINE, version=PQ_VERSION).compute()
    midpoint_df = pr.PyRanges(
        fragment_df.reset_index()[["chrom", "start", "end", "fragment_id"]]
        .assign(start=lambda x: ((x.start + x.end) * 0.5).round(0).astype(int))
        .eval("end = start + 1")
        .rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    fragment_to_bin = midpoint_df.join(bins, how="left").df[["fragment_id", "bin_id"]]
    fragment_to_bin = fragment_to_bin.set_index("fragment_id").sort_index()  # .astype(np.uint32)
    nulls = fragment_to_bin["bin_id"] == -1
    if nulls.any():
        logger.warning(
            "Some fragments did not overlap bins, removing from analysis:\n{}".format(
                fragment_to_bin[nulls].join(fragment_df)
            )
        )
        fragment_to_bin = fragment_to_bin[~nulls]
    binned_contacts = (
        contact_df.merge(fragment_to_bin, how="inner", right_index=True, left_on="align1_fragment_id")
        .merge(fragment_to_bin, how="inner", right_index=True, left_on="align2_fragment_id", suffixes=[None, "_2"])
        .rename(columns={"bin_id": "bin1_id", "bin_id_2": "bin2_id"})
    )
    pixels = binned_contacts.groupby(["bin1_id", "bin2_id"]).size().rename("count").astype(np.int32).reset_index()
    create_cooler(cooler_path, bins_df, pixels, ordered=True, symmetric_upper=True, ensure_sorted=True)
    c = Cooler(cooler_path)
    logger.info(f"Created cooler: {c.info}")
    return cooler_path


def export_to_paired_end_fastq(
    contact_table, output_prefix, reference_fasta, query, query_columns=None, read_length=50
):
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []

    columns.extend(
        [
            # "read_idx",
            "align1_fragment_id",
            "align1_chrom",
            "align1_fragment_start",
            "align1_fragment_end",
            "align1_strand",
            "align2_fragment_id",
            "align2_chrom",
            "align2_fragment_start",
            "align2_fragment_end",
            "align2_strand",
        ]
    )
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns)
    if query:
        query += " & (align1_fragment_id != align2_fragment_id) "
    else:
        query += "(align1_fragment_id != align2_fragment_id) "

    contact_df = contact_df.query(query)
    fastq1 = output_prefix + ".1.fastq"
    fastq2 = output_prefix + ".2.fastq"

    ref = FastaFile(reference_fasta)

    fhs = [open(fastq1, "w"), open(fastq2, "w")]
    qual_str = "I" * read_length
    contact_counter = defaultdict(int)
    for partition in range(contact_df.npartitions):
        df = (
            contact_df.get_partition(partition)
            .compute()
            .assign(
                pos1=lambda x: np.rint(x.eval("0.5 * (align1_fragment_start + align1_fragment_end)")).astype(int),
                pos2=lambda x: np.rint(x.eval("0.5 * (align2_fragment_start + align2_fragment_end)")).astype(int),
            )
        )
        for idx, row in df.iterrows():

            contact_counter[idx] += 1
            contact_idx = contact_counter[idx]
            read_name = f"read{idx:09}:{contact_idx:04}"

            read1 = ref.fetch(row.align1_chrom, row.pos1 - read_length, row.pos1)
            if row.align1_strand is False:
                read1 = revcomp(read1)
            fhs[0].write(f"@{read_name} 1:N:0:1\n{read1}\n{qual_str}\n")
            read2 = ref.fetch(row.align2_chrom, row.pos2 - read_length, row.pos2)
            if row.align2_strand is False:
                read2 = revcomp(read2)
            fhs[1].write(f"@{read_name} 2:N:0:1\n{read2}\n{qual_str}\n")
    return fastq1, fastq2
