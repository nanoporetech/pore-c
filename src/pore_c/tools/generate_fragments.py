from pysam import FastaFile
from typing import Pattern, List, NamedTuple, Iterator
import re

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

class SeqMatchPositions(NamedTuple):
    """Holds the results of a digestion"""
    seq_name: str
    positions: List[int]
    seq_length: int


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


def find_site_positions(regex: Pattern, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    return [m.start() for m in regex.finditer(seq.upper())]


def fragment_generator(reference_fasta: str, restriction_pattern: str) -> Iterator[SeqMatchPositions]:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""
    regex = create_regex(restriction_pattern)
    fasta_file = FastaFile(reference_fasta)

    for seqid in fasta_file.references:
        seq = fasta_file.fetch(seqid)
        yield SeqMatchPositions(seqid, find_site_positions(regex, seq), len(seq))
