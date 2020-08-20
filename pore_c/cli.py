import logging
from pathlib import Path

import click
import dask.dataframe as dd
import pandas as pd

import pore_c.catalogs as catalogs

from .catalogs import (
    RawReadCatalog,
    ReferenceGenomeCatalog,
    VirtualDigestCatalog,
)
from .cli_utils import (
    ExportDependentOption,
    NaturalOrderGroup,
    command_line_json,
    expand_output_prefix,
    filename_matches_regex,
    pipeable_sam_input,
    pipeable_sam_output,
)
from .config import INPUT_REFGENOME_REGEX, PQ_ENGINE, PQ_VERSION
from .settings import setup_logging
from .utils import DaskExecEnv


logger = setup_logging()


@click.group(cls=NaturalOrderGroup)
@click.option("-v", "--verbosity", count=True, help="Increase level of logging information, eg. -vvv")
@click.option("--quiet", is_flag=True, default=False, help="Turn off all logging", show_default=True)
@click.option("--dask-num-workers", type=int, default=1, help="Number of dask workers")
@click.option("--dask-use-threads", is_flag=True, default=False, help="Use threads instead of processes")
@click.option("--dask-threads-per-worker", type=int, default=1, help="Number of threads per worker")
@click.option(
    "--dask-scheduler-port",
    type=int,
    default=8786,
    help="The port to use for the dask scheduler, set to 0 to use a random port",
)
@click.option(
    "--dask-dashboard-port",
    type=int,
    default=8787,
    help="The port to use for the dask dashboard, set to 0 to use a random port",
)
@click.option("--dask-disable-dashboard", is_flag=True, default=False, help="Disable the dask dashboard")
@click.pass_context
def cli(
    ctx,
    verbosity,
    quiet,
    dask_num_workers,
    dask_use_threads,
    dask_threads_per_worker,
    dask_scheduler_port,
    dask_dashboard_port,
    dask_disable_dashboard,
):
    """Pore-C tools

    A suite of tools designed to analyse Oxford Nanopore reads with multiway chromatin contacts.
    """
    if quiet:
        logger.setLevel(logging.CRITICAL)
    elif verbosity > 0:
        LOG_LEVELS = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
        offset = 2
        idx = min(len(LOG_LEVELS) - 1, offset + verbosity)
        logger.setLevel(LOG_LEVELS[idx])
    else:
        logger.setLevel(logging.INFO)
    logger.debug("Logger set up")
    ctx.meta["dask_env"] = DaskExecEnv(
        n_workers=dask_num_workers,
        threads_per_worker=dask_threads_per_worker,
        processes=not dask_use_threads,
        scheduler_port=dask_scheduler_port,
        dashboard_port=None if dask_disable_dashboard else dask_dashboard_port,
    )


@cli.group(cls=NaturalOrderGroup, short_help="Pre-process reference genome files.")
@click.pass_context
def refgenome(ctx):
    pass


@refgenome.command(short_help="Pre-process a reference genome")
@click.argument("reference_fasta", type=click.Path(exists=True), callback=filename_matches_regex(INPUT_REFGENOME_REGEX))
@click.argument("output_prefix", callback=expand_output_prefix(ReferenceGenomeCatalog))
@click.option("--genome-id", type=str, help="An ID for this genome assembly")
@click.pass_context
def prepare(ctx, reference_fasta, output_prefix, genome_id):
    """Pre-process a reference genome for use by pore-C tools.

    Prepare a reference genome or draft assembly use by pore-c tools.
    This tool creates the following files:
    \b
        <output_prefix>.fa - A decompressed fasta file
        <output_prefix>.chromsizes - A tab-separated list of chromosome lengths
        <output_prefix>.metadata.csv - A csv with metadata about each of
        <output_prefix>.catalog.yaml - An intake catalog of these files

    """
    from pore_c.datasources import IndexedFasta
    import pandas as pd
    import re
    import subprocess as sp
    import pysam
    from shutil import copyfile

    logger.info(f"Adding reference genome under prefix: {output_prefix}")
    file_paths = ctx.meta["file_paths"]
    reference_fasta = Path(reference_fasta)

    parts_m = re.compile(INPUT_REFGENOME_REGEX).match(reference_fasta.name)
    stem, _, compression = parts_m.groups()
    if not genome_id:
        genome_id = stem

    dest_fasta = file_paths["fasta"]
    if compression == ".gz":
        comd = f"gunzip -cd {reference_fasta} > {dest_fasta}"
        logger.debug(f"Decompressing source fasta: {comd}")
        try:
            sp.check_call(comd, shell=True)
        except Exception as exc:  # noqa: F841
            logger.exception(f"Error creating bgzipped reference: {dest_fasta}")
            raise
    else:
        logger.debug("Copying {reference_fasta} to {dest_fasta}")
        copyfile(reference_fasta, dest_fasta)
    logger.debug("Creating faidx file")
    pysam.faidx(str(dest_fasta))

    ref_source = IndexedFasta(dest_fasta)
    ref_source.discover()
    chrom_lengths = {c["chrom"]: c["length"] for c in ref_source.metadata["chroms"]}
    chrom_df = pd.DataFrame(ref_source.metadata["chroms"])[["chrom", "length"]]
    chrom_df.to_csv(file_paths["chrom_metadata"], index=False)
    chrom_df.to_csv(file_paths["chromsizes"], sep="\t", header=None, index=False)
    metadata = {"chrom_lengths": chrom_lengths, "genome_id": genome_id}
    rg_cat = ReferenceGenomeCatalog.create(file_paths, metadata, {})
    logger.info("Added reference genome: {}".format(str(rg_cat)))


@refgenome.command(short_help="Virtual digest of a reference genome.")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("cut_on")
@click.argument("output_prefix", callback=expand_output_prefix(VirtualDigestCatalog))
@click.option(
    "--digest-type",
    type=click.Choice(["enzyme", "regex", "bin"]),
    default="enzyme",
    help="The type of digest you want to do",
)
@click.option("-n", "--n_workers", help="The number of dask_workers to use", default=1, show_default=True)
@click.pass_context
def virtual_digest(ctx, fasta, cut_on, output_prefix, digest_type, n_workers):
    """
    Carry out a virtual digestion of the genome/assembly FASTA.

    The DIGEST_TYPE sets what type of value you can use for CUT_ON.

    \b
        ---------------------------------------------------------------------
        | digest_type |  cut_on           | notes                           |
        ---------------------------------------------------------------------
        | enzyme      | NlaIII            | Enzyme name is case sensitive   |
        | enzyme      | HindIII           |                                 |
        | regex       | (GAATTC|GCGGCCGC) | Two site EcoRI and NotI         |
        | regex       | RAATY             | Degenerate site ApoI            |
        | bin         | 50000             | Create equal-width bins of 50k  |
        | bin         | 50k               | Create equal-width bins of 50k  |
        =====================================================================

    This tool will create the following output files:

    \b
        <output_prefix>.fragments.parquet

          A table containing the coordinates and some metadata on each fragment

        <output_prefix>.digest_stats.csv

          Per chromosome/contig summary statistics

        <output_prefix>.catalog.yaml

          An intake catalog of the digest files

    """
    from pore_c.analyses.reference import create_virtual_digest

    path_kwds = {key: val for key, val in ctx.meta["file_paths"].items() if key != "catalog"}
    logger.info(f"Creating virtual digest of {fasta}")
    frag_df, summary_stats = create_virtual_digest(fasta, digest_type, cut_on, **path_kwds)

    metadata = {"digest_type": digest_type, "digest_param": cut_on, "num_fragments": len(frag_df)}
    vd_cat = VirtualDigestCatalog.create(ctx.meta["file_paths"], metadata, {})
    logger.debug("Created Virtual Digest catalog: {}".format(vd_cat))


@refgenome.command(short_help="Create a hicRef file for a virtual digest.")
@click.argument("fragments_parquet", type=click.Path(exists=True))
@click.argument("hicref", type=click.Path(exists=False))
def fragments_to_hicref(fragments_parquet, hicref):
    """
    Convert a  *.fragments.parquet file to hicRef format
    """

    frag_df = pd.read_parquet(fragments_parquet, engine=PQ_ENGINE)
    with open(hicref, "w") as fh:
        for chrom, endpoints in frag_df.groupby("chrom")["end"].agg(lambda x: " ".join(map(str, x))).items():
            fh.write(f"{chrom} {endpoints}\n")

    logger.debug(f"Wrote hicRef file to {hicref}")


@cli.group(cls=NaturalOrderGroup, short_help="Analyse raw reads")
@click.pass_context
def reads(ctx):
    pass


@reads.command(short_help="Create a catalog file for a set of reads")  # noqa: F811
@click.argument("fastq", type=click.Path(exists=True))
@click.argument("output_prefix", callback=expand_output_prefix(RawReadCatalog))
@click.option(
    "--batch-size",
    help="The reads will be split into batches of this size for downstream processing",
    default=10_000,
    show_default=True,
)
@click.option("--min-read-length", help="The minimum length read to run through pore_c", default=1, show_default=True)
@click.option(
    "--max-read-length",
    help="The maximum length read to run through pore_c. Note that bwa mem can crash on very long reads",
    default=150_000,
    show_default=True,
)
@click.option("--min-qscore", help="The minimum read qscore", default=0, show_default=True)
@click.option("--max-qscore", help="The maximum read qscore", default=266, show_default=True)
@click.option("--user-metadata", callback=command_line_json, help="Additional user metadata to associate with this run")
@click.pass_context
def prepare(
    ctx, fastq, output_prefix, batch_size, min_read_length, max_read_length, min_qscore, max_qscore, user_metadata
):
    """Preprocess a set of reads for use with pore_c tools.

    This tool creates the following files:

    \b
        <output_prefix>.batch[batch_idx].fq.gz - Fastq with all the reads that pass the qscore and
          length filters. Fastqs are split so there are at most --batch-size reads per fastq
        <output_prefix>.fail.fq.gz - Reads that fail the filters
        <output_prefix>.read_metadata.parquet - Length and qscore metadata for each read and
          whether they pass the filter
        <output_prefix>.summary.csv - Summary stats for all/pass/fail reads
        <output_prefix>.catalog.yaml - An intake catalog

    """
    from pore_c.analyses.reads import prepare_fastq

    file_paths = ctx.meta["file_paths"]

    file_paths = catalogs.RawReadCatalog.generate_paths(output_prefix)
    path_kwds = {key: val for key, val in file_paths.items() if key != "catalog"}
    path_kwds["pass_fastq"] = Path(str(path_kwds["pass_fastq"]).replace(".batch.", ".batch{:d}."))
    summary = prepare_fastq(
        input_fastq=fastq,
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        min_qscore=min_qscore,
        max_qscore=max_qscore,
        chunksize=batch_size,
        **path_kwds,
    )

    catalog = RawReadCatalog.create(file_paths, {"summary_stats": summary}, user_metadata)
    logger.info("Created catalog for results: {}".format(catalog))


@cli.group(cls=NaturalOrderGroup, short_help="Analyse aligned pore_c reads")
def alignments():
    pass


@alignments.command(short_help="Reformat a BAM file to have a unique read name per alignment")
@click.argument("input_sam", type=str, callback=pipeable_sam_input)
@click.argument("output_sam", type=str, callback=pipeable_sam_output)
@click.option(
    "--input-is-bam",
    is_flag=True,
    default=False,
    is_eager=True,
    help="If piping a BAM from stdin (rather than sam)",
    show_default=True,
)
@click.option(
    "--output-is-bam",
    is_flag=True,
    default=False,
    is_eager=True,
    help="If piping a BAM to stdout (rather than sam)",
    show_default=True,
)
@click.option("--set-bx-flag", is_flag=True, default=False, help="Set the BX tag to the read name", show_default=True)
def reformat_bam(input_sam, output_sam, input_is_bam, output_is_bam, set_bx_flag):
    """Reformat query_name in INPUT_SAM  and write to OUTPUT_SAM

    This tool reformats an alignment file so that it works with downstream
    steps in the Pore-C pipeline. For both files you can supply '-' if you want
    to read/write from/to stdin/stdout. The 'query_name' field of the alignment
    file will be reformatted so that each alignment in the SAM file has a
    unique query name:

    \b
        <read_id> -> <read_id>:<read_idx>:<align_idx>

    Where 'read_idx' is a unique integer id for each read within the file and
    'align_idx' is a unique integer id for each alignment within the file. The
    tool also adds a 'BX' tag consisting of the 'read_id' to each record.

    """
    logger.debug(f"Reformatting alignments from {input_sam.filename} to {output_sam.filename}")
    read_indices = {}
    for align_idx, align in enumerate(input_sam.fetch(until_eof=True)):
        read_id = align.query_name
        read_idx = read_indices.get(read_id, None)
        if read_idx is None:
            read_idx = len(read_indices)
            read_indices[read_id] = read_idx
        if set_bx_flag:
            align.set_tag(tag="BX", value=align.query_name, value_type="Z")
        align.query_name = f"{read_id}:{read_idx}:{align_idx}"
        output_sam.write(align)
    output_sam.close()
    align_idx += 1
    num_reads = len(read_indices)
    logger.info(f"Processed {align_idx} alignments from {num_reads} reads")


@alignments.command(short_help="Parse a namesortd bam to pore-C alignment format")
@click.argument("input_bam", type=click.Path(exists=True))
@click.argument("output_table", type=click.Path(exists=False))
@click.option(
    "--alignment-haplotypes", type=click.Path(exists=True), help="The alignment to haplotype mapping from whatshap"
)
@click.pass_context
def create_table(ctx, input_bam, output_table, alignment_haplotypes):
    """Convert a BAM file to a tabular format sorted by read for downstream analysis

    """
    from pore_c import model
    from pysam import AlignmentFile

    tmp_table = output_table + ".tmp"
    logger.debug(f"Writing temporary unsorted data to {tmp_table}")
    af = AlignmentFile(input_bam)
    chrom_order = list(af.references)
    assert "NULL" not in chrom_order
    chrom_order.append("NULL")
    logger.debug(f"Chromosome order {chrom_order}")

    align_df = model.AlignmentRecord.to_dataframe(
        [model.AlignmentRecord.from_aligned_segment(a) for a in af], chrom_order=chrom_order
    )
    align_df = align_df.sort_values(["read_name"])
    num_aligns, num_reads = len(align_df), align_df.read_idx.nunique()
    logger.debug(f"Writing {num_aligns} alignments for {num_reads} reads to {output_table}")

    if alignment_haplotypes:
        ht_df = pd.read_csv(alignment_haplotypes, sep="\t")
        align_df = model.AlignmentRecord.update_dataframe_with_haplotypes(align_df, ht_df)

    align_df.to_parquet(output_table, engine=PQ_ENGINE, index=False, version=PQ_VERSION)
    g = align_df.groupby(["align_type"])
    summary = pd.concat({"num_reads": g["read_idx"].nunique(), "num_aligns": g.size()}).unstack(level=0)
    logger.info(f"Mapping summary:\n {summary}")

    haplotype_counts = (
        align_df.haplotype.value_counts()
        .rename_axis("haplotype")
        .to_frame()
        .rename(columns={"haplotype": "num_aligns"})
    )
    logger.info(f"Haplotype counts:\n {haplotype_counts}")


@alignments.command(short_help="Parse a namesortd bam to pore-C alignment format")
@click.argument("align_table", type=click.Path(exists=True))
@click.argument("fragments_table", type=click.Path(exists=True))
@click.argument("pore_c_table")
@click.option(
    "--mapping_quality_cutoff",
    type=int,
    default=1,
    help="Minimum mapping quality for an alignment to be considered",
    show_default=True,
)
@click.option(
    "--min_overlap_length",
    type=int,
    default=10,
    show_default=True,
    help="Minimum overlap in base pairs between an alignment and restriction fragment",
)
@click.option(
    "--containment_cutoff",
    type=float,
    default=99.0,
    show_default=True,
    help=(
        "Minimum percentage of a fragment included in an overlap for that "
        "fragment to be considered 'contained' within an alignment"
    ),
)
@click.pass_context
def assign_fragments(
    ctx, align_table, fragments_table, pore_c_table, mapping_quality_cutoff, min_overlap_length, containment_cutoff
):
    """For each alignment in ALIGN_TABLE either filter out or assign a fragment from FRAGMENT_TABLE


    """
    from pore_c.analyses.alignments import assign_fragments

    logger.info(f"Assigning fragments from {fragments_table} to alignments from {align_table}")
    align_table = pd.read_parquet(align_table, engine=PQ_ENGINE)
    fragment_df = pd.read_parquet(fragments_table, engine=PQ_ENGINE)
    pore_c_df = assign_fragments(
        align_table,
        fragment_df,
        mapping_quality_cutoff=mapping_quality_cutoff,
        min_overlap_length=min_overlap_length,
        containment_cutoff=containment_cutoff,
    )
    num_overlapping_alignments = len(pore_c_df)
    logger.info(f"Found {num_overlapping_alignments} overlapping alignments")
    pore_c_df.to_parquet(pore_c_table, engine=PQ_ENGINE, index=False, version=PQ_VERSION)
    logger.info(f"Fragments written to {pore_c_table}")
    logger.info("Summary: \n{}".format(pore_c_df.filter_reason.value_counts()))


@alignments.command(short_help="Parse a namesortd bam to pore-C alignment format")
@click.argument("pore_c_table", type=click.Path(exists=True))
@click.argument("output_pore_c_table", type=click.Path(exists=False))
@click.option(
    "--threshold",
    type=float,
    default=0.8,
    help="major:minor haplotype fraction must be greater than this value to assign a consensus",
    show_default=True,
)
@click.pass_context
def assign_consensus_haplotype(ctx, pore_c_table, output_pore_c_table, threshold):
    """Calculate a per-read consensus haplotype for each phase_set in ALIGN_TABLE and write the results back
    to OUTPUT_ALIGN_TABLE


    """
    from pore_c.model import PoreCRecord

    pore_c_table = pd.read_parquet(pore_c_table, engine=PQ_ENGINE)

    def get_most_freq_haplotype(_df):
        props = _df.div(_df.sum(axis=1), axis=0)
        res = pd.concat({"haplotype": props.idxmax(axis=1), "proportion": props.max(axis=1)}, axis=1)
        return res

    def _assign_consensus_haplotype(pore_c_table, threshold):
        chrom_dtype = pore_c_table.chrom.head(1).dtype
        meta = PoreCRecord.pandas_dtype(overrides={"chrom": chrom_dtype})

        pore_c_table = pore_c_table.set_index(["read_name", "chrom", "phase_set", "align_idx"])
        consensus_haplotypes = (
            pore_c_table.query("pass_filter == True and haplotype != -1")
            .groupby(level=["read_name", "chrom", "phase_set"], as_index=True)["haplotype"]
            .value_counts()  # count haplotypes
            .unstack(fill_value=0)  # convert to wide dataframe with one column per haplotype
            .pipe(get_most_freq_haplotype)  # get the most frequent haplotype by phase_set
            .query(f"({threshold} < proportion < 1.0)")  # if prop is 1 then we don't need to update
            .loc[:, ["haplotype"]]  # select the phase sets where consensus is possible
        )

        if len(consensus_haplotypes) > 0:
            # broadcast the consensus haplotype to all alignments for that read, phase_set combination
            consensus, _ = consensus_haplotypes.align(pore_c_table, join="left")
            # changes = (consensus != original).sum()
            # overwrite the haplotype values where consensus is possible
            pore_c_table.update(consensus, overwrite=True)
        pore_c_table = (
            pore_c_table.reset_index().astype(meta)[list(meta.keys())].sort_values(["read_name", "align_idx"])
        )
        return pore_c_table

    pore_c_df = _assign_consensus_haplotype(pore_c_table, threshold)

    pore_c_df.to_parquet(output_pore_c_table, engine=PQ_ENGINE, index=False, version=PQ_VERSION)


@alignments.command(short_help="Parses the alignment table and converts to pairwise contacts")
@click.argument("pore_c_table", type=click.Path(exists=True))
@click.argument("contact_table", type=click.Path(exists=False))
@click.pass_context
def to_contacts(ctx, pore_c_table, contact_table):
    """Covert the alignment table to a pairwise contact table and associated concatemer table

    """
    from pore_c.analyses.alignments import to_contacts

    pore_c_df = pd.read_parquet(pore_c_table, engine=PQ_ENGINE)

    contacts_df = to_contacts(pore_c_df).sort_values(["read_name"])
    num_contacts = len(contacts_df)
    contacts_df.to_parquet(contact_table, engine=PQ_ENGINE, index=False, version=PQ_VERSION)
    logger.info(f"Wrote {num_contacts} contacts to {contact_table}")
    if num_contacts == 0:
        logger.warning("No contacts found in {pore_c_table}")


@cli.group(cls=NaturalOrderGroup, short_help="Work the the contacts table")
def contacts():
    pass


@contacts.command(short_help="Summarise a contact table")
@click.argument("src_contact_tables", nargs=-1, type=click.Path(exists=True))
@click.argument("dest_contact_table", type=click.Path(exists=False))
@click.pass_context
def merge(ctx, src_contact_tables, dest_contact_table):
    import pyarrow.parquet as pq
    from pyarrow import dataset

    parts = []
    for i in src_contact_tables:
        md = pq.read_metadata(i)
        if md.num_rows == 0:
            logger.warning(f"The following contact file has no entries, removing from merge: {i}")
            continue
        parts.append(i)

    ds = dataset.dataset(parts, format="parquet")
    df = dd.read_parquet(parts, engine=PQ_ENGINE, version=PQ_VERSION, index=False)
    df.to_parquet(dest_contact_table, engine=PQ_ENGINE, version=PQ_VERSION, schema=ds.schema, write_index=False)


@contacts.command(short_help="Summarise a contact table")
@click.argument("contact_table", type=click.Path(exists=True))
@click.argument("read_summary_table", type=click.Path(exists=True))
@click.argument("concatemer_table", type=click.Path(exists=False))
@click.argument("concatemer_summary_csv", type=click.Path(exists=False))
@click.pass_context
def summarize(ctx, contact_table, read_summary_table, concatemer_table, concatemer_summary_csv):
    from pore_c.analyses.contacts import gather_concatemer_stats, summarize_concatemer_table
    from .model import PoreCConcatemerRecord

    concatemer_meta = PoreCConcatemerRecord.pandas_dtype()

    with ctx.meta["dask_env"]:
        contacts_df = dd.read_parquet(contact_table, engine=PQ_ENGINE, version=PQ_VERSION, index=False)
        concatemer_df = contacts_df.map_partitions(gather_concatemer_stats, meta=concatemer_meta).compute()
        concatemer_df.to_parquet(concatemer_table, engine=PQ_ENGINE, version=PQ_VERSION, index=False)

    long_summary_df = summarize_concatemer_table(concatemer_df, read_summary_table)
    long_summary_df.to_csv(concatemer_summary_csv)
    logger.info("Concatemer summary written to {}:\n {}".format(concatemer_summary_csv, long_summary_df.to_string()))


@contacts.command(short_help="Export contacts to various formats")
@click.argument("contact_table", type=click.Path(exists=True))
#  @click.argument("concatemer_table", type=click.Path(exists=True))
@click.argument("format", type=click.Choice(["cooler", "salsa_bed", "paired_end_fastq", "pairs"]))
@click.argument("output_prefix")
@click.option(
    "--min-mapping-quality",
    type=int,
    default=0,
    show_default=True,
    help="Both alignments have mapping qualities greater than this",
)
@click.option(
    "--min-align-base-qscore",
    type=int,
    default=0,
    show_default=True,
    help="Both alignments have mean base qualities greater than this",
)
@click.option("--cooler-resolution", help="The bin width of the resulting matrix", default=1000, show_default=True)
@click.option(
    "--fragment-table",
    cls=ExportDependentOption,
    export_formats=["cooler"],
    help="The fragment table for the corresponding virtual digest",
)
@click.option(
    "--by-haplotype",
    is_flag=True,
    help="Create a cooler for each pair of haplotypes (eg 1-1, 1-2, 2-2,...). Only valid with 'cooler'",
)
@click.option(
    "--chromsizes",
    cls=ExportDependentOption,
    export_formats=["cooler", "pairs"],
    help="The chromsizes file for the corresponding genome",
)
@click.option(
    "--reference-fasta",
    cls=ExportDependentOption,
    export_formats=["paired_end_fastq"],
    help="The reference genome used to generate the contact table",
)
@click.pass_context
def export(
    ctx,
    contact_table,
    format,
    output_prefix,
    min_mapping_quality,
    min_align_base_qscore,
    cooler_resolution,
    fragment_table,
    by_haplotype,
    chromsizes,
    reference_fasta,
):
    """
    Export contacts to the following formats:

     - cooler: a sparse representation of a contact matrix
     - paired_end_fastq: for each contact create a pseudo pair-end read using the reference genome sequence

    """
    from pore_c.analyses.contacts import (
        export_to_cooler,
        export_to_paired_end_fastq,
        export_to_salsa_bed,
        export_to_pairs,
    )

    columns = []
    query = []
    if min_mapping_quality:
        columns.extend(["align1_mapping_quality", "align2_mapping_quality"])
        query.append(
            f"(align1_mapping_quality > {min_mapping_quality}) & (align2_mapping_quality > {min_mapping_quality}) "
        )
    if min_align_base_qscore:
        columns.extend(["align1_align_base_qscore", "align2_align_base_qscore"])
        query.append(
            f"(align1_align_base_qscore > {min_align_base_qscore}) & "
            f"(align2_align_base_qscore > {min_align_base_qscore}) "
        )

    query = " & ".join(query)

    if format == "cooler":
        cooler_paths = export_to_cooler(
            contact_table,
            output_prefix,
            cooler_resolution,
            fragment_table,
            chromsizes,
            query,
            query_columns=columns,
            by_haplotype=by_haplotype,
        )
        for cooler_path in cooler_paths:
            logger.info(f"Wrote cooler to {cooler_path}")
    elif format == "salsa_bed":
        with ctx.meta["dask_env"]:
            bed_path = export_to_salsa_bed(contact_table, output_prefix, query, query_columns=columns)
        logger.info(f"Wrote cooler to {bed_path}")
    elif format == "paired_end_fastq":
        fastq1, fastq2 = export_to_paired_end_fastq(
            contact_table, output_prefix, reference_fasta, query, query_columns=columns
        )
        logger.info(f"Wrote reads to {fastq1} and {fastq2}")
    elif format == "pairs":
        pairs = export_to_pairs(contact_table, output_prefix, chromsizes, query, query_columns=columns)
        logger.info(f"Wrote contacts to {pairs}")
    else:
        raise NotImplementedError(format)


@cli.group(cls=NaturalOrderGroup, short_help="Misc tools")
def utils():
    pass


@utils.command()
@click.argument("input_parquet", type=click.Path(exists=True))
@click.argument("output_csv", type=click.Path(exists=False))
@click.pass_context
def parquet_to_csv(ctx, input_parquet, output_csv):
    """Convert a parquet file to CSV

    """

    with ctx.meta["dask_env"]:
        df = dd.read_parquet(input_parquet, engine=PQ_ENGINE, version=PQ_VERSION)
        df.to_csv(
            output_csv, single_file=True, compute=True, compression="gzip" if output_csv.endswith(".gz") else None
        )
