import click
import os.path
from pore_c.tools.generate_fragments import fragment_generator


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.
    """
    def list_commands(self, ctx):
        """List command names as they are in commands dict.
        """
        return self.commands.keys()


@click.group(cls=NaturalOrderGroup)
def cli():
    """Pore-C tools"""
    pass



@cli.command(short_help="Virtual digest of reference genome.")
@click.argument('reference_fasta', type=click.Path(exists=True))
@click.argument('restriction_pattern')
@click.argument('outfile', type=click.File('w'))
def generate_fragments(reference_fasta, restriction_pattern, outfile):
    """
    Virtual digest of RESTRICTION_PATTERN on REFERENCE_FASTA writing results to OUTFILE (use '-' for stdout).

    This is redundant with
    (https://github.com/aidenlab/juicer/blob/master/misc/generate_site_positions.py)
    in the juicer universe and generates an identical file. The included
    version of the tools supports degenerate sequences and sequence groups
    coded as REGEXes and prints to stdout.

    Some sample RESTRICTION_PATTERNs:

    \b
      - a single site (HindIII): "AAGCTT"
      - two sites (ecoRI and NotI): "(GAATTC|GCGGCCGC)"
      - degenerate site (BglI): "GCCNNNNNGGC"
      - degenerate site (ApoI): "RAATY"

    The OUTFILE is written with the following format:

    \b
    chr site1 site2 site3 ... chrLength
    \b
    as in:
    1 11160 12411 12461 ... 249250621
    2 11514 11874 12160 ... 243199373
    3 60138 60662 60788 ... 198022430

    """

    faidx_file = reference_fasta + ".fai"
    if not os.path.exists(faidx_file):
        raise IOError("Faidx file doesn't exist, please run 'samtools faidx {}'".format(reference_fasta))
    for seq_match_pos in fragment_generator(reference_fasta, restriction_pattern):
        outfile.write("{} {} {}\n".format(
            seq_match_pos.seq_name,
            " ".join(map(str, seq_match_pos.positions)),
            seq_match_pos.seq_length
        ))

@cli.command(short_help="Map the reads against a reference genome")
@click.argument('bam', type=click.Path(exists=True))
def map_reads(bam):
    print(bam)


@cli.command(short_help="Cluster mappings by read")
@click.argument('bam', type=click.Path(exists=True))
def cluster_reads(bam):
    print(bam)


@cli.command()
@click.argument('bam', type=click.Path(exists=True))
def map_to_fragments(bam):
    print(bam)

@cli.command()
@click.argument('bam', type=click.Path(exists=True))
def flatten_multiway(bam):
    print(bam)


@cli.command()
@click.argument('bam', type=click.Path(exists=True))
def hic_to_hicpro(bam):
    print(bam)

