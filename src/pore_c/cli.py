import logging
import os.path
from pathlib import Path

import click
import click_log
import pandas as pd
from intake import open_catalog

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.
        """
        return self.commands.keys()


@click.group(cls=NaturalOrderGroup)
@click_log.simple_verbosity_option(logger)
def cli():
    """Pore-C tools"""
    pass


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
    from intake.catalog.local import YAMLFileCatalog
    import pandas as pd

    logger.info("Adding reference genome under prefix: {}".format(output_prefix))
    fasta = Path(str(reference_fasta))
    try:
        stem, fasta_ext, compression_ext = fasta.name.split(".", 2)
        if not genome_id:
            genome_id = stem
    except:
        raise ValueError("Fasta file should be gzip compressed and should be in form {file_stem}.(fa|fasta|fna).gz")
    faidx_file = (fasta.parent / stem).with_suffix(".{}.{}.fai".format(fasta_ext, compression_ext))
    if not faidx_file.exists():
        raise IOError("Faidx file doesn't exist, please run 'samtools faidx {}'".format(reference_fasta))


    file_paths, exists = ReferenceGenomeCatalog.generate_paths(output_prefix)
    if exists:
        for file_id in exists:
            logger.error("Output file already exists for {}: {}".format(file_id, file_paths[file_id]))
        raise IOError()

    ref_source = IndexedFasta(fasta)
    ref_source.discover()
    chrom_lengths = {c['chrom']: c['length'] for c in ref_source.metadata['chroms']}
    chrom_df = pd.DataFrame(ref_source.metadata["chroms"])[["chrom", "length"]]
    chrom_df.to_csv(file_paths['chrom_metadata'], index=False)

    rg_cat = ReferenceGenomeCatalog.create(file_paths['catalog'], fasta, file_paths['chrom_metadata'], chrom_lengths, genome_id)
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
    digest_type, digest_param = cut_on.split(':')
    assert(digest_type in ['bin', 'enzyme', 'regex'])

    file_paths, exists = VirtualDigestCatalog.generate_paths(output_prefix)
    if exists:
        for file_id in exists:
            logger.error("Output file already exists for {}: {}".format(file_id, file_paths[file_id]))
        raise IOError()

    frag_df = create_virtual_digest(
        rg_cat.fasta,
        digest_type,
        digest_param,
        file_paths['fragments'],
        file_paths['digest_stats'],
        n_workers=n_workers
    )
    vd_cat = VirtualDigestCatalog.create(file_paths, Path(reference_catalog), digest_type, digest_param, len(frag_df))
    logger.debug("Created Virtual Digest catalog: {}".format(vd_cat))

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
    from pore_c.catalogs import AlignmentDfCatalog, ReferenceGenomeCatalog, VirtualDigestCatalog
    from pore_c.analyses.alignments import parse_alignment_bam

    file_paths, exists = AlignmentDfCatalog.generate_paths(output_prefix)
    if exists:
        for file_id in exists:
            logger.error("Output file already exists for {}: {}".format(file_id, file_paths[file_id]))
        raise IOError()

    vd_cat = open_catalog(str(virtual_digest_catalog))
    chrom_order = list(vd_cat.refgenome.metadata['chrom_lengths'].keys())
    fragment_df = vd_cat.fragments.read()
    final_stats = parse_alignment_bam(
        input_bam,
        fragment_df,
        alignment_table=file_paths['alignment'],
        read_table=file_paths['read'],
        overlap_table=file_paths['overlap'],
        alignment_summary=file_paths['alignment_summary'],
        read_summary=file_paths['read_summary'],
        n_workers=n_workers,
        chunksize=chunksize,
    )
    adf_cat = AlignmentDfCatalog.create(file_paths, Path(input_bam), Path(virtual_digest_catalog), final_stats)
    logger.info(str(adf_cat))

@cli.group(cls=NaturalOrderGroup, short_help="Convert between file formats")
def convert():
    pass


@convert.command(help="Convert from an alignment table to pairs format")
@click.argument("align_catalog", type=click.Path(exists=True))
@click.argument("reference_catalog", type=click.Path(exists=True))
@click.argument("pairs_file", type=click.Path(exists=False))
@click.option("-n", "--n_workers", help="The number of dask_workers to use", default=1)
def align_table_to_pairs(align_catalog, reference_catalog, pairs_file, n_workers):
    from pore_c.analyses.convert import convert_align_df_to_pairs

    align_cat = open_catalog(align_catalog)
    ref_cat = open_catalog(reference_catalog)
    ref_metadata = ref_cat.chrom_metadata.read()
    chrom_lengths = dict(ref_metadata[['chrom', 'length']].values)
    align_df = align_cat.align_table.to_dask()
    convert_align_df_to_pairs(align_df, chrom_lengths, Path(reference_catalog).name, pairs_file, n_workers=n_workers)


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
    #from pore_c.dashboard import main
    from bokeh.__main__ import main as bokeh_entry_point

    main_path = pore_c.__file__.rsplit("/", 1)[0] + "/dashboard/"
    sys.argv = ['bokeh', 'serve'] + [main_path] + list(bokeh_serve_args) + ["--args", align_catalog]
    bokeh_entry_point()


    #from bokeh.server.server import Server
    #from bokeh.command.util import build_single_handler_applications, die, report_server_init_errors
    #from tornado.ioloop import IOLoop

    #from pore_c.dashboard.main  import modify_doc
    #with report_server_init_errors():
    #    server = Server(
    #        {'/': modify_doc},
    #        io_loop = IOLoop(),
    #        address = "localhost",
    #        port = 8788
    #    )
    #    #server.start()
    #    server.run_until_shutdown()
    #    print(server)


@cli.command(short_help="create a .poreC file from a namesorted alignment of poreC data")
@click.argument("filter_catalog", type=click.Path(exists=True))
@click.argument("fragment_map_catalog", type=click.Path(exists=True))
def map_to_fragments(filter_catalog, fragment_map_catalog):
    from intake import open_catalog
    from pore_c.datasources import IndexedBedFile
    from ncls import NCLS
    import numpy as np

    filter_cat = open_catalog(filter_catalog)
    fragment_cat = open_catalog(fragment_map_catalog)

    fragment_df = fragment_cat.fragment_df.read().set_index(["fragment_id"])  # .sort_values(['chrom', 'start'])
    fragment_intervals = {}
    for chrom, chrom_df in fragment_df.groupby("chrom"):
        if chrom == "NULL":
            continue
        fragment_intervals[chrom] = NCLS(chrom_df.start.values, chrom_df.end.values, chrom_df.index.values)

    for chunk in filter_cat.align_table.read_chunked():
        for chrom, chrom_df in chunk.groupby("chrom"):
            if chrom in fragment_intervals:
                chunk_indices, frag_indices = fragment_intervals[chrom].all_overlaps_both(
                    chrom_df.start.astype(int).values,
                    chrom_df.end.astype(int).values,
                    chrom_df.index.astype(int).values,
                )
                if len(chunk_indices) == 0:
                    continue
                df_a = chrom_df.reindex(chunk_indices).reset_index().rename(columns={"index": "hit_idx"})
                df_b = (
                    fragment_df.reindex(frag_indices)
                    .reset_index()
                    .rename(columns={"index": "fragment_id", "start": "fragment_start", "end": "fragment_end"})
                    .loc[:, ["fragment_start", "fragment_end", "fragment_id"]]
                )
                overlap_df = (
                    df_a.join(df_b)
                    .assign(
                        overlap_start=lambda x: np.where(x.fragment_start > x.start, x.fragment_start, x.start).astype(
                            int
                        ),
                        overlap_end=lambda x: np.where(x.fragment_end < x.end, x.fragment_end, x.end).astype(int),
                    )
                    .eval("overlap_length = overlap_end - overlap_start")
                )

                print(overlap_df.head())
        break
    raise ValueError(filter_cat)
    map_to_fragments_tool(input_bed, fragment_bed_file, output_porec, method, stats)
