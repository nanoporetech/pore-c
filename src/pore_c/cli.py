import click
import click_log
import logging
import os.path
from pore_c.tools.generate_fragments import create_fragment_map
from pore_c.tools.generate_fragments import create_bin_file as create_bin_file_tool
from pore_c.tools.map_to_bins import bin_hic_data as bin_hic_data_tool
from pore_c.tools.cluster_reads import cluster_reads as cluster_reads_tool
from pore_c.tools.cluster_reads import measure_overlaps as measure_overlaps_tool
from pore_c.tools.cluster_reads import remove_contained_segments as remove_contained_segments_tool
from pore_c.tools.map_to_frags import map_to_fragments as map_to_fragments_tool
from pore_c.tools.poreC_flatten import flatten_multiway as flatten_multiway_tool

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


@cli.command(short_help="Virtual digest of reference genome.")
@click.argument('reference_fasta', type=click.Path(exists=True))
@click.argument('restriction_pattern')
@click.argument('bedfile', type=click.Path(exists=False))
@click.argument('hicref', type=click.Path(exists=False))
def generate_fragments(reference_fasta, restriction_pattern, bedfile, hicref):
    """
    Carry out a virtual digestion of RESTRICTION_PATTERN on REFERENCE_FASTA writing results to BEDFILE and HICREF.

    The RESTRICTION_PATTERN can be specified using either the name of a restriction enzyme as available in the
    Biopython Restriction package (eg. HindIII, case sensitive) or as a python regular expression prefixed
    by 'regex:'. Note that the positions returned by these two methods will differ as the biopython method
    will return fragments delimited by the actual cut site, whereas the regular expression will return
    the positions of the recognition patterns.

    Some sample RESTRICTION_PATTERNs:

    \b
      - an enzyme name, note case sensitive: "HindIII"
      - a single site (HindIII): "regex:AAGCTT"
      - two sites (ecoRI and NotI): "regex:(GAATTC|GCGGCCGC)"
      - degenerate site (BglI): "regex:GCCNNNNNGGC"
      - degenerate site (ApoI): "regex:RAATY"

    The BEDFILE contains an entry for each of the fragments generated by the by the digest.

    \b
    chr\tstart\tend\tfragment_id
    \b

    The HICREF follows the Juicer format generated by this script:
    (https://github.com/aidenlab/juicer/blob/master/misc/generate_site_positions.py):

    \b
    chr site1 site2 site3 ... chrLength
    \b
    as in:
    1 11160 12411 12461 ... 249250621
    2 11514 11874 12160 ... 243199373
    3 60138 60662 60788 ... 198022430

    """
    logger.info("Creating fragments from reference genome: {} and digestion pattern {}".format(reference_fasta, restriction_pattern))
    faidx_file = reference_fasta + ".fai"
    if not os.path.exists(faidx_file):
        raise IOError("Faidx file doesn't exist, please run 'samtools faidx {}'".format(reference_fasta))
    frag_map = create_fragment_map(reference_fasta, restriction_pattern)
    logger.debug("Created FragmentMap: {}".format(frag_map))
    frag_map.save_to_bed(bedfile)
    logger.info("FragmentMap saved to bed file: {}".format(bedfile))
    frag_map.save_to_HiCRef(hicref)
    logger.info("FragmentMap saved to HicREF file: {}".format(hicref))

@cli.command(short_help="Cluster mappings by read")
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('keep_bam', type=click.Path(exists=False))
@click.argument('discard_bam', type=click.Path(exists=False))
@click.option("--trim", default=20, type=int, help="The number of bases to ignore at the end of each aligned segment when calculating overlaps.")
@click.option("--mapping_quality_cutoff", default=0, type=int, help="The minimum mapping quality for a alignment segment to be kept")
@click.option('--alignment_stats', default = None, type=click.Path(exists=False), help="A filename for storing logged data about fragment assignment on a per-alignment basis.")
@click.option("--contained", is_flag=True, help="If the contained flag is raised, cluster-filtering is not done, but instead alignments are removed based on whether they are fully contained in another alignment.") 
def cluster_reads(input_bam, keep_bam, discard_bam, trim, contained, mapping_quality_cutoff, alignment_stats):
    if contained:
        if alignment_stats is not None:
            remove_contained_segments_tool(input_bam, keep_bam, discard_bam, mapping_quality_cutoff, alignment_stats)
        else:
            remove_contained_segments_tool(input_bam, keep_bam, discard_bam, mapping_quality_cutoff)
    else:
        if alignment_stats is not None:
            num_reads, num_reads_kept, num_aligns, num_aligns_kept = cluster_reads_tool(input_bam, keep_bam, discard_bam, trim, mapping_quality_cutoff, alignment_stats)
        else:
            num_reads, num_reads_kept, num_aligns, num_aligns_kept = cluster_reads_tool(input_bam, keep_bam, discard_bam, trim, mapping_quality_cutoff)


@cli.command(short_help="create a .poreC file from a namesorted alignment of poreC data")
@click.argument('input_bed', type=click.Path(exists=True))# bedtools overlap file between the filtered bam and the reference fragment file
@click.argument('fragment_bed_file',type=click.Path(exists=True))# A reference fragment bed file generated by the pore_c generate-fragments command.
@click.argument('output_porec', type=click.Path(exists=False))
@click.option('--method', default = 'overlap', type = str, help="The method of determining fragment mapping of an alignment")
@click.option("--stats", default=None, type=click.Path(exists=False), help="A filename for storing the per-mapping logging data about fragment mapping.")
def map_to_fragments(input_bed, fragment_bed_file, output_porec, method, stats):
    map_to_fragments_tool(input_bed,fragment_bed_file, output_porec, method, stats)


@cli.command(short_help="In each sequencing read, this tool measures the overlap intervals between all pairs of alignments (for diagnosing alignment filtering).") 
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_table', type=click.Path(exists=False))
@click.option('--no_zero', is_flag = True, help="for pairs of alignments that do not overlap, do not include a table entry. This cuts down the table size dramatically.")
def measure_overlaps(input_bam, output_table, no_zero):
    measure_overlaps_tool(input_bam,output_table, no_zero)

@cli.command(short_help = "Flatten down a pore-C file filled with multiway contacts to a single specified contact dimension." )
@click.argument('input_porec',type=click.Path(exists=True))
@click.argument('output_porec', type=click.Path(exists=False))
@click.option('--sort', is_flag = True, help="Sort the monomers within each contact according to fragment ID. This does not sort the entries as they might need to be sorted for conversion to other formats or for visualisation tools..")
@click.option('--size', default=2, type=int, help="The size of the generated contacts in the output file. Default is 2.")
def flatten_multiway(input_porec, output_porec,size,sort):
    flatten_multiway_tool(input_porec, output_porec,size,sort)

@cli.command(short_help = "Flatten down a pore-C file filled with multiway contacts to a single specified contact dimension." )
@click.argument('input_fai',type=click.Path(exists=True))
@click.argument('output_bin_bed', type=click.Path(exists=False))
@click.argument('size', default=1000000, type=int)# help="The bin size for the file. Default is 10**6 bp."
def create_bin_file(input_fai, output_bin_bed, size):
    create_bin_file_tool(input_fai, output_bin_bed,size)

#@click.argument('size', default=1000000, type=int, help="The bin size for the file. Default is 10**6 bp.") #this might be unnecessary as an argument, since it's implied by the reference bedfile.

@cli.command()
@click.argument('hictxt', type=click.Path(exists=True))
@click.argument('hic_ref', type=click.Path(exists=True))
@click.argument('bin_ref', type=click.Path(exists=True))
@click.argument('output_bin_matrix', type=click.Path(exists=False))
def bin_hic_data(hictxt,hic_ref,bin_ref,output_bin_matrix):
    bin_hic_data_tool(hictxt,output_bin_matrix, hic_ref, bin_ref)

