import click
import os.path
#from pore_c.tools.generate_fragments import fragment_generator
from pore_c.tools.generate_fragments import digest

from pore_c.tools.cluster_reads import cluster_reads as cluster_reads_tool
from pore_c.tools.map_to_frags import map_to_fragments as map_to_fragments_tool
from pore_c.tools.poreC_flatten import flatten_multiway as flatten_multiway_tool


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
def biogenerate_fragments(reference_fasta, restriction_pattern, outfile):
    #if restriction pattern contains multiple patterns,
    #   this line turns that into a RE list
    #   if it is a single pattern, it properly
    #   packages it up as a RE list used by the functions
    #   on up the stack
    re = restriction_pattern.split(' ')
    for entry in digest(re, reference_fasta):
        outfile.write(entry + "\n")


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


@cli.command(short_help="Cluster mappings by read")
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('keep_bam', type=click.Path(exists=False))
@click.argument('discard_bam', type=click.Path(exists=False))
@click.option("--trim", default=20, type=int, help="The number of bases to ignore at the end of each aligned segment when calculating overlaps")
@click.option("--mapping_quality_cutoff", default=0, type=int, help="The minimum mapping quality for a alignment segment to be kept")
def cluster_reads(input_bam, keep_bam, discard_bam, trim, mapping_quality_cutoff):
    num_reads, num_reads_kept, num_aligns, num_aligns_kept = cluster_reads_tool(input_bam, keep_bam, discard_bam, trim, mapping_quality_cutoff)
    print(num_reads, num_reads_kept, num_aligns, num_aligns_kept)


@cli.command(short_help="create a .poreC file from a namesorted alignment of poreC data")
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('reference',type=click.Path(exists=True))
@click.argument('output_porec', type=click.Path(exists=False))
@click.option('--method', default = 'start', type = str, help="The method of determining fragment mapping of an alignment")
def map_to_fragments(input_bam, reference, output_porec, method):
    map_to_fragments_tool(input_bam,reference, output_porec, method)


@cli.command(short_help = "Flatten down a pore-C file filled with multiway contacts to a single specified contact dimension." )
@click.argument('input_porec',type=click.Path(exists=True))
@click.argument('output_porec', type=click.Path(exists=False))
@click.option('--sort', default=False, type=bool, help="Sort the monomers in each contact according to fragment ID.")
@click.option('--size', default=2, type=int, help="The size of the generated contacts in the output file. Default is 2.")
def flatten_multiway(input_porec, output_porec,size,sort):
    flatten_multiway_tool(input_porec, output_porec,size,sort)

@cli.command()
@click.argument('bam', type=click.Path(exists=True))
def hic_to_hicpro(bam):
    print(bam)

