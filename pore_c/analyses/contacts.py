import logging
import os
import shutil
import subprocess as sp
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
    contact_table,
    output_prefix,
    cooler_resolution,
    fragment_table,
    chromsizes,
    query,
    query_columns=None,
    by_haplotype=False,
):

    results = []
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []
    columns.extend(["align1_fragment_id", "align2_fragment_id"])
    if by_haplotype:
        columns.extend(["align1_haplotype", "align2_haplotype"])
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns)
    if query:
        contact_df = contact_df.query(query)

    chrom_dict = pd.read_csv(
        chromsizes, sep="\t", header=None, names=["chrom", "size"], index_col=["chrom"], squeeze=True
    )
    # create even-widht bins using cooler
    bins_df = binnify(chrom_dict, cooler_resolution)
    bins_df.index.name = "bin_id"
    # convert to ranges for overlap
    bins = pr.PyRanges(bins_df.reset_index().rename(columns={"start": "Start", "end": "End", "chrom": "Chromosome"}))

    fragment_df = dd.read_parquet(fragment_table, engine=PQ_ENGINE, version=PQ_VERSION).compute()
    midpoint_df = pr.PyRanges(
        fragment_df.reset_index()[["chrom", "start", "end", "fragment_id"]]
        .assign(start=lambda x: ((x.start + x.end) * 0.5).round(0).astype(int))
        .eval("end = start + 1")
        .rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    # use a pyranges joing to assign fragments to bins
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

    # use a join to assign each end of a contact to a bin
    binned_contacts = (
        contact_df.merge(fragment_to_bin, how="inner", right_index=True, left_on="align1_fragment_id")
        .merge(fragment_to_bin, how="inner", right_index=True, left_on="align2_fragment_id", suffixes=[None, "_2"])
        .rename(columns={"bin_id": "bin1_id", "bin_id_2": "bin2_id"})
    )

    if not by_haplotype:
        cooler_path = output_prefix + ".cool"
        # group size == number of contacts per bin_pair
        pixels = binned_contacts.groupby(["bin1_id", "bin2_id"]).size().rename("count").astype(np.int32).reset_index()
        create_cooler(cooler_path, bins_df, pixels, ordered=True, symmetric_upper=True, ensure_sorted=True)
        c = Cooler(cooler_path)
        logger.info(f"Created cooler: {c.info}")
        results.append(cooler_path)
    else:
        tmp_parquet = output_prefix + ".tmp.pq"
        pixels = (
            # create a key to groupy by haplotype pair, order of haplotypes doesn't matter
            binned_contacts.assign(
                hap_key=lambda x: x[["align1_haplotype", "align2_haplotype"]].apply(
                    lambda y: "{}_{}".format(*sorted(y)).replace("-1", "nohap"), axis=1, meta="object"
                )
            )
            .groupby(["hap_key", "bin1_id", "bin2_id"])
            .size()
            .rename("count")
            .astype(np.int32)
            .reset_index()
            .astype({"hap_key": "category"})
        )

        # save to a temporary parquet file, this might not be necessary
        # but want to avoid the whole contact matrix hitting memory
        pixels.to_parquet(
            tmp_parquet,
            write_metadata_file=True,
            partition_on=["hap_key"],
            write_index=False,
            engine=PQ_ENGINE,
            version=PQ_VERSION,
        )

        pixels = dd.read_parquet(tmp_parquet, engine=PQ_ENGINE, version=PQ_VERSION, columns=["hap_key"])
        hap_keys = pixels["hap_key"].unique().compute()
        # create a cooler for each haplotype pair
        for hap_key in hap_keys:
            cooler_path = f"{output_prefix}.{hap_key}.cool"
            pixels = dd.read_parquet(
                tmp_parquet,
                filters=[("hap_key", "==", hap_key)],
                index=False,
                engine=PQ_ENGINE,
                version=PQ_VERSION,
                columns=["bin1_id", "bin2_id", "count"],
            )
            create_cooler(cooler_path, bins_df, pixels, ordered=True, symmetric_upper=True, ensure_sorted=True)
            c = Cooler(cooler_path)
            logger.info(f"Created cooler: {c.info}")
            results.append(cooler_path)

        shutil.rmtree(tmp_parquet)

    return results


def export_to_salsa_bed(contact_table, output_prefix, query, query_columns):
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []

    columns.extend(
        [
            "align1_chrom",
            "align1_start",
            "align1_end",
            "align1_strand",
            "align1_mapping_quality",
            "align2_chrom",
            "align2_start",
            "align2_end",
            "align2_strand",
            "align2_mapping_quality",
        ]
    )
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns)
    if query:
        contact_df = contact_df.query(query)
    bed_file = output_prefix + ".salsa.bed"

    dtypes = contact_df.head(1).dtypes
    strand_dtype = pd.CategoricalDtype(["+", "-"], ordered=False)
    meta = {
        "chrom": dtypes["align1_chrom"],
        "start": dtypes["align1_chrom"],
        "end": dtypes["align1_chrom"],
        "read_pair_id": str,
        "strand": strand_dtype,
        "mapping_quality": dtypes["align1_mapping_quality"],
    }

    def to_salsa_long(df):
        contact_idx = df.groupby(level="read_idx")["align1_chrom"].transform(lambda x: np.arange(len(x), dtype=int))
        df["align1_contact_idx"] = contact_idx
        df["align2_contact_idx"] = contact_idx
        new_columns = pd.MultiIndex.from_tuples(
            [c.replace("align", "").split("_", 1) for c in df.columns], names=["pair_idx", "value"]
        )
        df.columns = new_columns
        df = df.stack(level="pair_idx").reset_index()
        df["read_pair_id"] = df[["read_idx", "contact_idx", "pair_idx"]].apply(
            lambda x: f"read{x.read_idx:012}_{x.contact_idx}/{x.pair_idx}", axis=1
        )
        df["strand"] = df["strand"].replace({True: "+", False: "-"}).astype(strand_dtype)
        return df[["chrom", "start", "end", "read_pair_id", "strand", "mapping_quality"]]

    contact_df.map_partitions(to_salsa_long, meta=meta).to_csv(
        bed_file, single_file=True, sep="\t", header=False, index=False
    )

    return bed_file


def export_to_pairs(contact_table, output_prefix, chromsizes, query, query_columns):
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []

    columns.extend(
        [
            "contact_is_direct",
            "contact_read_distance",
            "align1_align_idx",
            "align1_chrom",
            "align1_fragment_start",
            "align1_fragment_end",
            "align1_strand",
            "align1_mapping_quality",
            "align2_align_idx",
            "align2_chrom",
            "align2_fragment_start",
            "align2_fragment_end",
            "align2_strand",
            "align2_mapping_quality",
        ]
    )
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns)
    chrom_dict = pd.read_csv(
        chromsizes, sep="\t", header=None, names=["chrom", "size"], index_col=["chrom"], squeeze=True
    )

    if query:
        contact_df = contact_df.query(query)

    dtypes = contact_df.head(1).dtypes
    strand_dtype = pd.CategoricalDtype(["+", "-"], ordered=False)
    junction_dtype = pd.CategoricalDtype(["DJ", "IJ"], ordered=False)
    meta = {
        "readID": str,
        "chr1": dtypes["align1_chrom"],
        "pos1": dtypes["align1_fragment_start"],
        "chr2": dtypes["align2_chrom"],
        "pos2": dtypes["align2_fragment_start"],
        "strand1": strand_dtype,
        "strand2": strand_dtype,
        "pair_type": junction_dtype,
        "align1_idx": dtypes["align1_align_idx"],
        "align2_idx": dtypes["align2_align_idx"],
        "distance_on_read": dtypes["contact_read_distance"],
    }

    pairs_file = output_prefix + ".pairs"
    # dask to_csv doesn't support append writing, so need to concatenate header after the fact
    header_file = output_prefix + ".pairs.header"
    body_file = output_prefix + ".pairs.body"

    header_lines = lines = ["## pairs format 1.0", "#shape: upper triangle"]
    header_lines.extend(["#chromsize: {} {}".format(*_) for _ in chrom_dict.items()])
    header_lines.append("#columns: {}".format(" ".join(meta.keys())))

    with open(header_file, "w") as output_fh:
        output_fh.write("{}\n".format("\n".join(lines)))

    def to_pairs_df(df):
        df["contact_idx"] = df.groupby(level="read_idx")["align1_chrom"].transform(
            lambda x: np.arange(len(x), dtype=int)
        )
        df = (
            df.reset_index()
            .assign(
                readID=lambda x: x[["read_idx", "contact_idx"]]
                .astype(int)
                .apply(lambda y: f"read{y.read_idx:012}_{y.contact_idx:d}", axis=1),
                pos1=lambda x: np.rint(x.eval("0.5 * (align1_fragment_start + align1_fragment_end)")).astype(int),
                pos2=lambda x: np.rint(x.eval("0.5 * (align2_fragment_start + align2_fragment_end)")).astype(int),
                align1_strand=lambda x: x["align1_strand"].replace({True: "+", False: "-"}).astype(strand_dtype),
                align2_strand=lambda x: x["align2_strand"].replace({True: "+", False: "-"}).astype(strand_dtype),
                pair_type=lambda x: x["contact_is_direct"].replace({True: "DJ", False: "IJ"}).astype(junction_dtype),
            )
            .rename(
                columns=dict(
                    align1_chrom="chr1",
                    align2_chrom="chr2",
                    align1_strand="strand1",
                    align2_strand="strand2",
                    align2_align_idx="align2_idx",
                    align1_align_idx="align1_idx",
                    contact_read_distance="distance_on_read",
                )
            )
        )
        return df[list(meta.keys())]

    contact_df.map_partitions(to_pairs_df, meta=meta).to_csv(
        body_file, mode="w+", single_file=True, sep="\t", header=False, index=False
    )

    sp.check_call(f"cat {header_file} {body_file} > {pairs_file}", shell=True)
    os.unlink(header_file)
    os.unlink(body_file)

    return pairs_file


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
