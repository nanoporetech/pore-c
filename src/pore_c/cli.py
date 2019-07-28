import json
import logging
from pathlib import Path

import click
from intake import open_catalog

from .settings import setup_logging
import  pore_c.catalogs as catalogs

logger = setup_logging()
logger.warning("here")


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.
        """
        return self.commands.keys()


def command_line_json(ctx, param, value):
    # TODO: add support for json from file
    try:
        res = json.loads(value)
    except Exception as exc:  # noqa: F841
        logger.exception("Not valid json")
        raise
    return res


@click.group(cls=NaturalOrderGroup)
@click.option("-v", "--verbosity", count=True, help="Increase level of logging information, eg. -vvv")
@click.option("--quiet", is_flag=True, help="Turn off all logging")
def cli(verbosity, quiet):
    """Pore-C tools"""
    if quiet:
        logger.setLevel(logging.CRITICAL)
    elif verbosity > 0:
        LOG_LEVELS = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
        offset = 2
        idx = min(len(LOG_LEVELS) - 1, offset + verbosity)
        logger.setLevel(LOG_LEVELS[idx])
    else:
        pass
    logger.debug("Logger set up")


@cli.group(cls=NaturalOrderGroup, short_help="Operations on the reference genome sequence")
def refgenome():
    pass


@refgenome.command(short_help="Create a catalog file for the reference genome.")
@click.argument("reference_fasta", type=click.Path(exists=True))
@click.argument("output_prefix")
@click.option("--genome-id", type=str, help="An ID for this genome assembly")
def catalog(reference_fasta, output_prefix, genome_id=None):
    """Preprocess a reference genome for use by pore_c tools
    """
    from pore_c.datasources import IndexedFasta
    from pore_c.catalogs import ReferenceGenomeCatalog
    import pandas as pd

    logger.info("Adding reference genome under prefix: {}".format(output_prefix))
    fasta = Path(str(reference_fasta))
    try:
        stem, fasta_ext, compression_ext = fasta.name.split(".", 2)
        if not genome_id:
            genome_id = stem
    except Exception as e:
        raise ValueError(
            "Fasta file should be gzip compressed and should be in form {file_stem}.(fa|fasta|fna).gz\n{}".format(e)
        )
    faidx_file = (fasta.parent / stem).with_suffix(".{}.{}.fai".format(fasta_ext, compression_ext))
    if not faidx_file.exists():
        raise IOError("Faidx file doesn't exist, please run 'samtools faidx {}'".format(reference_fasta))

    file_paths = ReferenceGenomeCatalog.generate_paths(output_prefix)
    path_kwds = {key: val for key, val in file_paths.items() if key != 'catalog'}

    ref_source = IndexedFasta(fasta)
    ref_source.discover()
    chrom_lengths = {c["chrom"]: c["length"] for c in ref_source.metadata["chroms"]}
    chrom_df = pd.DataFrame(ref_source.metadata["chroms"])[["chrom", "length"]]
    chrom_df.to_csv(file_paths["chrom_metadata"], index=False)
    chrom_df.to_csv(file_paths["chromsizes"], sep="\t", header=None, index=False)
    metadata = {'chrom_lengths': chrom_lengths, 'genome_id': genome_id}
    file_paths['fasta'] = fasta
    rg_cat = ReferenceGenomeCatalog.create(file_paths, metadata, {})
    logger.info("Added reference genome: {}".format(str(rg_cat)))


@refgenome.command(short_help="Virtual digest of reference genome.")
@click.argument("reference_catalog", type=click.Path(exists=True))
@click.argument("cut_on")
@click.argument("output_prefix")
@click.option("-n", "--n_workers", help="The number of dask_workers to use", default=1)
def virtual_digest(reference_catalog, cut_on, output_prefix, n_workers):
    """
    Carry out a virtual digestion of CUT_ON on REFERENCE_FASTA writing results to BEDFILE and HICREF.

    \b

    Some sample CUT_ONs:

    \b
      - an enzyme name, note case sensitive: "enzyme:HindIII"
      - a fixed width bin: "bin:50k"
      - a single site (HindIII): "regex:AAGCTT"
      - two sites (ecoRI and NotI): "regex:(GAATTC|GCGGCCGC)"
      - degenerate site (BglI): "regex:GCCNNNNNGGC"
      - degenerate site (ApoI): "regex:RAATY"

    """
    from pore_c.catalogs import ReferenceGenomeCatalog, VirtualDigestCatalog
    from pore_c.analyses.reference import create_virtual_digest

    rg_cat = ReferenceGenomeCatalog(reference_catalog)
    digest_type, digest_param = cut_on.split(":")
    assert digest_type in ["bin", "enzyme", "regex"]

    file_paths = VirtualDigestCatalog.generate_paths(output_prefix)
    path_kwds = {key: val for key, val in file_paths.items() if key != 'catalog'}

    frag_df = create_virtual_digest(
        rg_cat.fasta,
        digest_type,
        digest_param,
        n_workers=n_workers,
        **path_kwds,
    )

    metadata = {
        'digest_type': digest_type,
        'digest_param': digest_param,
        'num_fragments': len(frag_df),
    }
    file_paths['refgenome_catalog'] = Path(reference_catalog)
    vd_cat = VirtualDigestCatalog.create(file_paths, metadata, {})
    logger.debug("Created Virtual Digest catalog: {}".format(vd_cat))


@cli.group(cls=NaturalOrderGroup, short_help="Analyse raw reads")
def reads():
    pass


@reads.command(short_help="Create a catalog file for a set of reads")  # noqa: F811
@click.argument("fastq", type=click.Path(exists=True))
@click.argument("output_prefix")
@click.option("--min-read-length", help="The minimum length read to run through porec", default=1)
@click.option(
    "--max-read-length",
    help="The maximum length read to run through porec. Note that bwa mem can crash on very long reads",
    default=500000,
)
@click.option("--user-metadata", callback=command_line_json, help="Additional user metadata to associate with this run")
def catalog(fastq, output_prefix, min_read_length, max_read_length, user_metadata):
    """Preprocess a reference genome for use by pore_c tools
    """
    from pore_c.analyses.reads import filter_fastq

    file_paths = catalogs.RawReadCatalog.generate_paths(output_prefix)
    path_kwds = {key: val for key, val in file_paths.items() if key != 'catalog'}
    summary = filter_fastq(
        input_fastq=fastq,
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        **path_kwds
    )

    catalog = catalogs.RawReadCatalog.create(file_paths, {'summary_stats': summary}, user_metadata)
    logger.info("Created catalog for results: {}".format(catalog))

    c1 = open_catalog(str(file_paths['catalog']))
    logger.info(c1)

@cli.group(cls=NaturalOrderGroup, short_help="Analyse aligned porec reads")
def alignments():
    pass


@alignments.command(short_help="Parse a namesortd bam to pore-C alignment format")
@click.argument("input_bam", type=click.Path(exists=True))
@click.argument("virtual_digest_catalog", type=click.Path(exists=True))
@click.argument("output_prefix")
@click.option("-n", "--n_workers", help="The number of dask_workers to use", default=1)
@click.option("--chunksize", help="Number of reads per processing chunk", default=50000)
def parse(input_bam, virtual_digest_catalog, output_prefix, n_workers, chunksize):
    """Filter the read-sorted alignments in INPUT_BAM and save the results under OUTPUT_PREFIX

    """
    from pore_c.catalogs import AlignmentDfCatalog
    from pore_c.analyses.alignments import parse_alignment_bam

    file_paths, exists = AlignmentDfCatalog.generate_paths(output_prefix)
    if exists:
        for file_id in exists:
            logger.error("Output file already exists for {}: {}".format(file_id, file_paths[file_id]))
        raise IOError()

    vd_cat = open_catalog(str(virtual_digest_catalog))
    fragment_df = vd_cat.fragments.read()
    final_stats = parse_alignment_bam(
        input_bam,
        fragment_df,
        alignment_table=file_paths["alignment"],
        read_table=file_paths["read"],
        overlap_table=file_paths["overlap"],
        alignment_summary=file_paths["alignment_summary"],
        read_summary=file_paths["read_summary"],
        n_workers=n_workers,
        chunksize=chunksize,
    )
    adf_cat = AlignmentDfCatalog.create(file_paths, Path(input_bam), Path(virtual_digest_catalog), final_stats)
    logger.info(str(adf_cat))


@cli.group(cls=NaturalOrderGroup, short_help="Create pairs files")
def pairs():
    pass


@pairs.command(help="Convert from an alignment table to pairs format")
@click.argument("align_catalog", type=click.Path(exists=True))
@click.argument("output_prefix")
@click.option("-n", "--n_workers", help="The number of dask_workers to use", default=1)
def from_alignment_table(align_catalog, output_prefix, n_workers):
    from pore_c.analyses.pairs import convert_align_df_to_pairs
    from pore_c.catalogs import PairsFileCatalog

    file_paths, exists = PairsFileCatalog.generate_paths(output_prefix)
    if exists:
        for file_id in exists:
            logger.error("Output file already exists for {}: {}".format(file_id, file_paths[file_id]))
        raise IOError()

    adf_cat = open_catalog(align_catalog)
    rg_cat = adf_cat.virtual_digest.refgenome
    chrom_lengths = rg_cat.metadata["chrom_lengths"]
    genome_id = rg_cat.metadata["genome_id"]
    align_df = adf_cat.alignment.to_dask()
    convert_align_df_to_pairs(align_df, chrom_lengths, genome_id, file_paths["pairs"], n_workers=n_workers)

    pair_cat = PairsFileCatalog.create(file_paths, Path(align_catalog))
    logger.info(str(pair_cat))


@cli.group(cls=NaturalOrderGroup, short_help="Dashboard")
def dashboard():
    pass


@dashboard.command(help="Alignment dashboard", context_settings=dict(ignore_unknown_options=True))
@click.argument("align_catalog", nargs=1, type=click.Path(exists=True))
@click.argument("bokeh_serve_args", nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def alignment(ctx, align_catalog, bokeh_serve_args):
    import sys
    import pore_c

    # from pore_c.dashboard import main
    from bokeh.__main__ import main as bokeh_entry_point

    main_path = pore_c.__file__.rsplit("/", 1)[0] + "/dashboard/"
    sys.argv = ["bokeh", "serve"] + [main_path] + list(bokeh_serve_args) + ["--args", align_catalog]
    bokeh_entry_point()
