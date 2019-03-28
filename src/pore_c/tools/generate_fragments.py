from pysam import FastaFile
from typing import Pattern, List, NamedTuple, Iterator
import re

from pore_c.model import SeqDigest, FragmentMap

# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    'ACGTWSMKRYBDHVNacgtwsmkrybdhvn-)(][',
    'TGCAWSKMYRVHDBNtgcawskmyrvhdbn-()[]'
)

# translate degenerate bases to regex set
DEGENERATE_TRANS = str.maketrans(
    dict([
        ("N","[ACGT]"),
        ("V","[ACG]"),
        ("H","[ACT]"),
        ("D","[AGT]"),
        ("B","[CGT]"),
        ("W","[AT]"),
        ("S","[CG]"),
        ("M","[AC]"),
        ("K","[GT]"),
        ("R","[AG]"),
        ("Y","[CT]"),
    ])
)



def revcomp(seq: str) -> str:
    """Return the reverse complement of a string:
    """
    return seq[::-1].translate(COMPLEMENT_TRANS)

def replace_degenerate(pattern: str) -> str:
    """Replace degenerate bases with regex set"""
    return pattern.translate(DEGENERATE_TRANS)


def create_regex(pattern: str) -> Pattern :
    """Takes a raw restriction digest site in the form of a regular expression string and returns a
    regular expression object consisting of both forward and reverse complement versions of the pattern"""
    ###

    site_raw = pattern.replace('(',' ').replace(')',' ').replace(' ','').replace('|',' ').split()
    sites_raw = []
    for entry in sites_raw:
        sites_raw.append(revcomp(entry))

    sites_raw += site_raw
    if len(sites_raw) > 1:
        fwd_rev_pattern = "(" + '|'.join(sorted(list(set(sites_raw)))) + ")"
    else:
        fwd_rev_pattern = "(" + sites_raw[0] + ")"

    ###
    fwd_rev_pattern = replace_degenerate(fwd_rev_pattern)

    try:
        regex = re.compile(fwd_rev_pattern)
    except:
        raise ValueError("Error compiling regex for pattern {}, redundance form: {}".format(pattern, fwd_rev_pattern))
    return regex



def find_site_positions(regex_or_enzyme: str, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    if regex_or_enzyme.startswith('regex:'):
        regex = create_regex(regex_or_enzyme.split('regex:', 1)[1])
        return find_site_positions_regex(regex, seq)
    else:
        return find_site_positions_biopython(regex_or_enzyme, seq)


def find_site_positions_regex(regex: Pattern, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    return [m.start() for m in regex.finditer(seq.upper())]


def find_site_positions_biopython(enzyme: str, seq:str) -> List[int]:
    from Bio import Restriction
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        raise ValueError("Enzyme not found: {}".format(enzyme))
    #print("{} cut site: {}".format(enz, enz.elucidate()))
    s = Seq(seq, IUPACAmbiguousDNA())
    positions = [_ -1 for _ in enz.search(s)]
    if len(positions) == 0:
        return [0]
    else:
        return positions


def fragment_generator(reference_fasta: str, restriction_pattern: str) -> Iterator[SeqDigest]:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""
    fasta_file = FastaFile(reference_fasta)

    for seqid in fasta_file.references:
        seq = fasta_file.fetch(seqid)
        yield SeqDigest(seqid, find_site_positions(restriction_pattern, seq), len(seq))


def create_fragment_map(reference_fasta: str, restriction_pattern: str) -> FragmentMap:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""
    frag_map = FragmentMap.from_digest_iter(
        fragment_generator(reference_fasta, restriction_pattern)
    )
    return frag_map


#a tool for creating bin level partitionings of the reference genome
def create_bin_file(reference_fai: str, bedfile_out: str,  bin_size: int) -> None:
    """generate a bedfile from a reference fai that describes the bin tiling that covers the entire genome."""
    frag_num = 0
    f_out = open(bedfile_out,'w')
    for line in open(reference_fai):
        l = line.strip().split()
        l[1] = int(l[1])
        frags = list(range(0, bin_size * (1+(l[1] // bin_size)), bin_size))
        frags = list(zip(frags[:-1],frags[1:]))
        tail_frag = (((l[1] -1 ) // bin_size) * bin_size , l[1])
        frags.append(tail_frag)
        for idx, entry in enumerate( frags ) :
            f_out.write('{ch}\t{st}\t{en}\t{num}\n'.format(ch = l[0] , st = entry[0], en = entry[1], num = idx + frag_num))
        frag_num += len(frags)
        


