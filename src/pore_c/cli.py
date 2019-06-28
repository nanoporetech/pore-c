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


@cli.command(short_help="Create a catalog file for the reference genome.")
@click.argument("reference_fasta", type=click.Path(exists=True))
@click.argument("output_prefix")
def add_refgenome(reference_fasta, output_prefix):
    """Preprocess a reference genome for use by pore_c tools
    """
    from pore_c.datasources import IndexedFasta
    from intake.catalog.local import YAMLFileCatalog
    import pandas as pd
    import yaml

    logger.info("Adding reference genome under prefix: {}".format(output_prefix))
    fasta = Path(str(reference_fasta))
    try:
        stem, fasta_ext, compression_ext = fasta.name.split(".", 2)
    except:
        raise ValueError(
            "Fasta file should be gzip compressed and should be in form {file_stem}.(fa|fasta|fna).gz"
        )

    faidx_file = (fasta.parent / stem).with_suffix(".{}.{}.fai".format(fasta_ext, compression_ext))
    if not faidx_file.exists():
        raise IOError(
            "Faidx file doesn't exist, please run 'samtools faidx {}'".format(reference_fasta)
        )
    output_prefix = Path(output_prefix)
    catalog_path = output_prefix.with_suffix(".refgenome.catalog.yaml")
    if catalog_path.exists():
        raise IOError("Catalog file exists: {}".format(catalog_path))

    ref_source = IndexedFasta(fasta)
    ref_source.discover()
    chrom_df = pd.DataFrame(ref_source.metadata["chroms"])[["chrom", "length"]].assign(
        organism="UNKNOWN", molecule_type="chromosome|plasmid"
    )
    metadata_csv = output_prefix.with_suffix(".metadata.csv")
    chrom_df.to_csv(metadata_csv, index=False)

    catalog_data = {
        "name": "pore_c_reference_genome",
        "description": "A reference genome for use with pore-c tools",
        "sources": {
            "fasta": {
                "driver": "pore_c.datasources.IndexedFasta",
                "args": {"urlpath": "{}".format(fasta.resolve())},
            },
            "chrom_metadata": {
                "driver": "csv",
                "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(metadata_csv.name)},
            },
        },
    }
    with catalog_path.open("w") as fh:
        fh.write(yaml.dump(catalog_data))
    cat = YAMLFileCatalog(str(catalog_path))


@cli.command(short_help="Virtual digest of reference genome.")
@click.argument("reference_catalog", type=click.Path(exists=True))
@click.argument("restriction_pattern")
@click.argument("output_prefix")
def virtual_digest(reference_catalog, restriction_pattern, output_prefix):
    """
    Carry out a virtual digestion of RESTRICTION_PATTERN on REFERENCE_FASTA writing results to BEDFILE and HICREF.

    \b

    Some sample RESTRICTION_PATTERNs:

    \b
      - an enzyme name, note case sensitive: "HindIII"
      - a fixed width bin: "bin:50k"
      - a single site (HindIII): "regex:AAGCTT"
      - two sites (ecoRI and NotI): "regex:(GAATTC|GCGGCCGC)"
      - degenerate site (BglI): "regex:GCCNNNNNGGC"
      - degenerate site (ApoI): "regex:RAATY"

    """
    from pore_c.analyses.reference import create_virtual_digest

    logger.info(
        "Creating fragments from reference genome: {} and digestion pattern {}".format(
            reference_catalog, restriction_pattern
        )
    )
    reference_cat = open_catalog(reference_catalog)
    output_prefix = Path(output_prefix)
    cat = create_virtual_digest(reference_cat.fasta, restriction_pattern, output_prefix)
    logger.debug("Created Virtual Digest catalog: {}".format(cat.path))


@cli.command(short_help="Filter alignments")
@click.argument("input_bam", type=click.Path(exists=True))
@click.argument("output_prefix")
@click.option(
    "--save_fail_bam", is_flag=True, help="Save alignments that fail filter to {prefix}.pass.bam"
)
@click.option(
    "--save_align_table",
    is_flag=True,
    help="Save the filter status of each alignment to {prefix}.alignment.parquet",
)
@click.option(
    "--save_read_table", is_flag=True, help="Save read-level statistics to {prefix}.read.parquet"
)
def filter_alignments(input_bam, output_prefix, save_fail_bam, save_align_table, save_read_table):
    from pore_c.analyses.alignments import filter_alignments as filt
    from pathlib import Path

    output_prefix = Path(output_prefix)

    pass_bam = output_prefix.with_suffix(".pass.bam")
    if not pass_bam.parent.exists():
        pass_bam.parent.mkdir()
    fail_bam = output_prefix.with_suffix(".fail.bam") if save_fail_bam else None
    align_table = output_prefix.with_suffix(".alignment.parquet") if save_align_table else None
    read_table = output_prefix.with_suffix(".read.parquet") if save_read_table else None
    catalog_file = output_prefix.with_suffix(".catalog.yaml")

    fail = False
    for outfile in [pass_bam, fail_bam, align_table, read_table]:
        if outfile is not None and outfile.exists():
            logger.error("Output file already exists: {}".format(outfile))
            fail = True

    if fail:
        raise click.ClickException("An error was encountered while setting up alignment filters")

    res = filt(
        input_bam,
        pass_bam=pass_bam,
        fail_bam=fail_bam,
        align_table=align_table,
        read_table=read_table,
        catalog_file=catalog_file,
    )
    logger.info(res)


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

    fragment_df = fragment_cat.fragment_df.read().set_index(
        ["fragment_id"]
    )  # .sort_values(['chrom', 'start'])
    fragment_intervals = {}
    for chrom, chrom_df in fragment_df.groupby("chrom"):
        if chrom == "NULL":
            continue
        fragment_intervals[chrom] = NCLS(
            chrom_df.start.values, chrom_df.end.values, chrom_df.index.values
        )

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
                df_a = (
                    chrom_df.reindex(chunk_indices)
                    .reset_index()
                    .rename(columns={"index": "hit_idx"})
                )
                df_b = (
                    fragment_df.reindex(frag_indices)
                    .reset_index()
                    .rename(
                        columns={
                            "index": "fragment_id",
                            "start": "fragment_start",
                            "end": "fragment_end",
                        }
                    )
                    .loc[:, ["fragment_start", "fragment_end", "fragment_id"]]
                )
                overlap_df = (
                    df_a.join(df_b)
                    .assign(
                        overlap_start=lambda x: np.where(
                            x.fragment_start > x.start, x.fragment_start, x.start
                        ).astype(int),
                        overlap_end=lambda x: np.where(
                            x.fragment_end < x.end, x.fragment_end, x.end
                        ).astype(int),
                    )
                    .eval("overlap_length = overlap_end - overlap_start")
                )

                print(overlap_df.head())
        break
    raise ValueError(filter_cat)
    map_to_fragments_tool(input_bed, fragment_bed_file, output_porec, method, stats)
