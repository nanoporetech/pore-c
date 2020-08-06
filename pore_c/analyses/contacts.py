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

from pore_c.analyses.reads import read_length_stats
from pore_c.analyses.reference import revcomp
from pore_c.config import PQ_ENGINE, PQ_VERSION, SHORT_RANGE_CUTOFF
from pore_c.model import (
    PoreCConcatemerRecord,
    PoreCConcatemerRecordDf,
    PoreCContactRecordDf,
)


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
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns, index=False)
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

        pixels = dd.read_parquet(tmp_parquet, engine=PQ_ENGINE, version=PQ_VERSION, columns=["hap_key"], index=False)
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
            "read_name",
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
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns, index=False)
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
        df = df.set_index("read_name")
        contact_idx = df.groupby(level="read_name")["align1_chrom"].transform(lambda x: np.arange(len(x), dtype=int))
        df["align1_contact_idx"] = contact_idx
        df["align2_contact_idx"] = contact_idx
        new_columns = pd.MultiIndex.from_tuples(
            [c.replace("align", "").split("_", 1) for c in df.columns], names=["pair_idx", "value"]
        )
        df.columns = new_columns
        df = df.stack(level="pair_idx").reset_index()
        df["read_pair_id"] = df[["read_name", "contact_idx", "pair_idx"]].apply(
            lambda x: f"read{x.read_name}_{x.contact_idx}/{x.pair_idx}", axis=1
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
            "read_name",
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
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns, index=False)
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
        df["contact_idx"] = df.groupby("read_name")["align1_chrom"].transform(lambda x: np.arange(len(x), dtype=int))
        df = (
            df.reset_index()
            .assign(
                readID=lambda x: x[["read_name", "contact_idx"]].apply(
                    lambda y: f"read{y.read_name}_{y.contact_idx:d}", axis=1
                ),
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
            "read_name",
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
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns, index=False)
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

        df["pos1"].clip(lower=read_length, inplace=True)
        df["pos2"].clip(lower=read_length, inplace=True)
        for idx, row in df.iterrows():

            contact_counter[row.read_name] += 1
            contact_idx = contact_counter[row.read_name]
            read_name = f"{row.read_name}:{contact_idx:04}"

            read1 = ref.fetch(row.align1_chrom, row.pos1 - read_length, row.pos1)
            if row.align1_strand is False:
                read1 = revcomp(read1)
            fhs[0].write(f"@{read_name} 1:N:0:1\n{read1}\n{qual_str}\n")
            read2 = ref.fetch(row.align2_chrom, row.pos2 - read_length, row.pos2)
            if row.align2_strand is False:
                read2 = revcomp(read2)
            fhs[1].write(f"@{read_name} 2:N:0:1\n{read2}\n{qual_str}\n")
    return fastq1, fastq2


def export_to_merged_no_dups(contact_table, output_prefix, reference_fasta, query, query_columns=None, read_length=50):
    if query_columns:
        columns = query_columns[:]
    else:
        columns = []

    columns.extend(
        [
            "read_name",
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
    contact_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, columns=columns, index=False)
    if query:
        query += " & (align1_fragment_id != align2_fragment_id) "
    else:
        query += "(align1_fragment_id != align2_fragment_id) "

    contact_df = contact_df.query(query)
    mnd_file = output_prefix + ".mnd.txt"
    fh = open(mnd_file, "w")

    ref = FastaFile(reference_fasta)
    cigar_str = f"{read_length}M"
    contact_counter = defaultdict(int)
    for partition in range(contact_df.npartitions):
        df = (
            contact_df.get_partition(partition)
            .compute()
            .astype({"align1_strand": "uint8", "align2_strand": "uint8"})
            .assign(
                pos1=lambda x: np.rint(x.eval("0.5 * (align1_fragment_start + align1_fragment_end)")).astype(int),
                pos2=lambda x: np.rint(x.eval("0.5 * (align2_fragment_start + align2_fragment_end)")).astype(int),
            )
            .drop_duplicates(subset=["pos1", "pos1"])
        )

        df["pos1"].clip(lower=read_length, inplace=True)
        df["pos2"].clip(lower=read_length, inplace=True)
        for _, row in df.iterrows():

            contact_counter[row.read_name] += 1
            contact_idx = contact_counter[row.read_name]
            read_name = f"{row.read_name}:{contact_idx:04}"

            read1 = ref.fetch(row.align1_chrom, row.pos1 - read_length, row.pos1)
            if row.align1_strand == 1:
                read1 = revcomp(read1)
            read2 = ref.fetch(row.align2_chrom, row.pos2 - read_length, row.pos2)
            if row.align2_strand == 1:
                read2 = revcomp(read2)
            # <str1> <chr1> <pos1> <frag1>
            # <str2> <chr2> <pos2> <frag2>
            # <mapq1> <cigar1> <sequence1>
            # <mapq2> <cigar2> <sequence2>
            # <readname1> <readname2>
            fh.write(
                f"{row.align1_strand} {row.align1_chrom} {row.pos1} {row.align1_fragment_id} "
                f"{row.align2_strand} {row.align2_chrom} {row.pos2} {row.align2_fragment_id} "
                f"60 {cigar_str} {read1} "
                f"60 {cigar_str} {read2} "
                f"{read_name}:1 {read_name}:2\n"
            )

    return mnd_file


def gather_concatemer_stats(contact_df: PoreCContactRecordDf) -> PoreCConcatemerRecordDf:

    contact_df["is_short_range"] = contact_df["contact_fragment_distance"] < SHORT_RANGE_CUTOFF
    by_read = contact_df.groupby("read_name", as_index=True)
    read_stats = by_read[["read_length", "read_idx"]].first()

    num_fragments = (
        by_read[["align1_fragment_id", "align2_fragment_id"]]
        .agg(lambda x: set(x.unique()))
        .apply(lambda x: len(x.align1_fragment_id.union(x.align2_fragment_id)), axis=1)
        .rename("num_fragments")
    )

    num_reads = len(by_read)
    by_read_and_type = contact_df.groupby(["read_name", "contact_is_direct"], as_index=True)

    def flatten_index(df, sep="_"):
        df.columns = [sep.join(t) for t in df.columns.to_flat_index()]
        return df

    contact_counts = (
        contact_df.groupby(["read_name", "contact_is_direct", "contact_is_cis"])
        .size()
        .rename(index={True: "direct", False: "indirect"}, level="contact_is_direct")
        .rename(index={True: "cis", False: "trans"}, level="contact_is_cis")
        .unstack(fill_value=0, level="contact_is_direct")
        .unstack(fill_value=0, level="contact_is_cis")
        .pipe(flatten_index)
        .astype(int)
        .eval("total_cis = direct_cis + indirect_cis")
        .eval("total_trans = direct_trans + indirect_trans")
        .eval("direct = direct_cis + direct_trans")
        .eval("indirect = indirect_cis + indirect_trans")
        .eval("total = total_cis + total_trans")
        .add_suffix("_contacts")
        .eval("read_order = direct_contacts + 1")
    )

    short_range_counts = (
        contact_df.loc[contact_df.contact_is_cis]
        .groupby(["read_name", "contact_is_direct", "is_short_range"], as_index=True)
        .size()
        .rename(index={True: "short_range", False: "long_range"}, level=-1)
        .rename(index={True: "direct", False: "indirect"}, level=-2)
        .unstack(fill_value=0, level=-2)
        .unstack(fill_value=0, level=-1)
        .pipe(flatten_index)
        .eval("total_short_range = direct_short_range + indirect_short_range")
        .eval("total_long_range = direct_long_range + indirect_long_range")
        .add_suffix("_cis_contacts")
    )

    haplotype_stats = by_read["haplotype_pair_type"].value_counts().unstack(fill_value=0)
    drop = []
    for cat in contact_df.haplotype_pair_type.cat.categories:
        if cat == "null" or cat == "trans":
            if cat in haplotype_stats.columns:
                drop.append(cat)
        elif cat not in haplotype_stats.columns:
            haplotype_stats[cat] = 0
    if drop:
        haplotype_stats = haplotype_stats.drop(columns=drop)
    haplotype_stats = haplotype_stats.add_prefix("haplotype_")

    max_distance = (
        by_read_and_type[["contact_genome_distance", "contact_fragment_distance"]]
        .max()
        .unstack(fill_value=0)
        .astype(int)
        .rename(columns={True: "direct", False: "indirect"})
    )
    max_distance.columns = ["max_{1}_{0}".format(*_) for _ in max_distance.columns]
    dtype = PoreCConcatemerRecord.pandas_dtype()
    res = (
        contact_counts.join(haplotype_stats, how="left")
        .join(read_stats, how="left")
        .join(max_distance, how="left")
        .join(short_range_counts, how="left")
        .join(num_fragments, how="left")
        .fillna({c: 0 for c in short_range_counts.columns})
        .reset_index()
        # .astype(dtype)
        [list(dtype.keys())]
    )
    assert len(res) == num_reads
    return res


def summarize_concatemer_table(concatemer_df: PoreCConcatemerRecordDf, read_summary_table: str) -> pd.DataFrame:
    concatemer_read_length_stats = read_length_stats(concatemer_df["read_length"])
    # read_subset_dtype = pd.CategoricalDtype(["all", "pass", "fail"])
    read_summary_table = (
        pd.read_csv(read_summary_table)
        .query("read_subset == 'all'")
        .assign(index=0)
        # .astype({"read_subset": read_subset_dtype})
        .set_index(["index", "read_subset"])
        .unstack(fill_value=0)
        .rename({"all": "all_reads", "pass": "pass_read_qc", "fail": "fail_read_qc"}, axis=1, level=1)
    )

    # print(read_summary_table)

    contact_counts = (
        concatemer_df.loc[:, lambda x: [c for c in x.columns if c.endswith("_contacts")]]
        .sum()
        .astype(int)
        .rename("count")
        .to_frame()
    )
    contact_counts.index = pd.MultiIndex.from_tuples(
        [c.split("_", 1) for c in contact_counts.index], names=["subset", "value"]
    )
    contact_props = (
        contact_counts.div(contact_counts.xs("contacts", level="value", drop_level=True), level="subset", axis=0)
        .mul(100.0)
        .rename(columns={"count": "perc"})
    )
    contact_section = (
        contact_counts.join(contact_props)
        .rename_axis("variable", axis=1)
        .stack()
        .rename("value")
        .to_frame()
        .reorder_levels([0, 2, 1])
        .sort_index(level=[0, 1], sort_remaining=False)
    )

    read_counts = read_summary_table.xs("num_sequences", axis=1).copy()
    read_data = {"all": read_counts.loc[0, "all_reads"], "concatemer": len(concatemer_df)}
    read_data["fail_filter"] = read_data["all"] - read_data["concatemer"]
    read_data = pd.Series(read_data, name="count").to_frame()
    read_data["perc"] = 100.0 * read_data["count"] / read_data.loc["all", "count"]

    read_data = (
        read_data.rename_axis("variable", axis=1).stack().rename("value").to_frame().reorder_levels([1, 0]).sort_index()
    )

    bases = read_summary_table.xs("total_bases", axis=1).copy()
    base_data = {"all": bases.loc[0, "all_reads"], "concatemer": concatemer_read_length_stats["total_bases"]}
    base_data["fail_filter"] = base_data["all"] - base_data["concatemer"]
    base_data = pd.Series(base_data, name="bases").to_frame()
    base_data["perc"] = 100.0 * base_data["bases"] / base_data.loc["all", "bases"]
    base_data = (
        base_data.rename_axis("variable", axis=1).stack().rename("value").to_frame().reorder_levels([1, 0]).sort_index()
    )
    read_length_data = {
        "all": read_summary_table.at[0, ("N50", "all_reads")],
        "concatemer": concatemer_read_length_stats["N50"],
    }
    read_length_data = (
        pd.Series(read_length_data, name="value")
        .to_frame()
        .rename_axis("subset")
        .assign(variable="N50")
        .set_index("variable", append=True)
        .reorder_levels([1, 0])
        .sort_index(level=0)
    )

    # fragments_per_alignment = concatemer_df.eval("(1.0 * num_fragments) /  read_order")
    # print(fragments_per_alignment.describe())
    # fragment_data = {
    #    "mean": fragments_per_alignment.mean(),
    #    "median": fragments_per_alignment.median(),
    # }
    # print(fragment_data)

    length_bins = pd.IntervalIndex.from_breaks([1, 2, 3, 4, 6, 11, 21, 50, int(1e9)])
    length_bin_labels = {}
    for i in length_bins:
        if i.length == 1:
            label = str(i.right)
        elif i.length >= int(1e8):
            label = f"gt_{i.left}"
        else:
            label = "{}-{}".format(i.left + 1, i.right)
        length_bin_labels[i] = label

    read_order_hist = (
        pd.cut(concatemer_df["read_order"], length_bins, labels=length_bin_labels)
        .value_counts()
        .rename("count")
        .sort_index()
    )
    read_order_hist = read_order_hist.to_frame().rename(index=length_bin_labels).rename_axis("concatemer_order")
    read_order_hist["perc"] = read_order_hist.div(read_order_hist["count"].sum(), axis=1) * 100
    read_order_hist = (
        read_order_hist.rename_axis("variable", axis=1)
        .stack()
        .rename("value")
        .to_frame()
        .reorder_levels([1, 0])
        .sort_index(level="variable", sort_remaining=False)
    )

    total_contacts = contact_section.at[("total", "count", "contacts"), "value"]
    density_data = (
        pd.Series(
            {
                "all": 1e9 * total_contacts / base_data.at[("bases", "all"), "value"],
                "concatemer": 1e9 * total_contacts / base_data.at[("bases", "concatemer"), "value"],
            }
        )
        .rename_axis("divisor")
        .rename("value")
        .to_frame()
        .assign(variable="contacts_per_gb")
        .set_index("variable", append=True)
        .reorder_levels([1, 0])
    )

    def normalize_levels(dfs):
        max_levels = max([len(df.index.names) for df in dfs.values()])
        level_labels = [f"level_{x}" for x in range(max_levels)]
        res = {}
        for key, val in dfs.items():
            existing_levels = val.index.names
            add_levels = max_levels - len(existing_levels)
            rename_levels = dict(zip(existing_levels, level_labels[add_levels:]))
            for x in range(add_levels):
                val[level_labels[x]] = ""
                val = val.set_index(level_labels[x], append=True)
            new_labels = [rename_levels.get(l, l) for l in val.index.names]
            val.index.rename(new_labels, inplace=True)
            val = val.reorder_levels(level_labels)
            res[key] = val
        res = pd.concat(res, axis=0, names=["section"])
        idx_names = res.index.names
        res = res.T.convert_dtypes().astype(str).T
        res.index.rename(idx_names, inplace=True)
        return res

    long_df = normalize_levels(
        {
            "reads": read_data,
            "read_length": read_length_data,
            "bases": base_data,
            "contacts": contact_section,
            "density": density_data,
            "concatemer_order": read_order_hist,
        }
    )
    return long_df
