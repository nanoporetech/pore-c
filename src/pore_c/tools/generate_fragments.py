#from pysam import FastaFile
from typing import Pattern, List, NamedTuple, Iterator
#import re

import gzip
from Bio import Seq,SeqIO
from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict

#takes in the sequence argument provided by name or by sequence,
#  looks through the Bio.Restriction class members
#  for a matching enzyme.
# returns a restriction batch with either a single enzyme or multiple 
# enzymes depending on what was specified
def loadRestrictionEnzymes(rest_dict, enzyme_targets):
    #generate a table of all restriction enzymes in Biopython
    RE_ref = {}
    for name, data in rest_dict.items():
        enzyme = getattr(Restriction, name)
        RE_ref[name] = enzyme
        RE_ref[data["site"]] = enzyme
    #identify enzymes by motif or by name
    rb = Restriction.RestrictionBatch()
    for arg in enzyme_targets:
        if arg in RE_ref.values():
            rb.add(getattr(Restriction,arg))
        elif arg in RE_ref.keys():
            rb.add(RE_ref[arg])
    return rb


def digest(enzyme_targets, ref):
    restriction_batch = loadRestrictionEnzymes(rest_dict,enzyme_targets)
#    raise TypeError(type(restriction_batch))
    if '.gz' in ref:
        fIn = gzip.open(ref,'rt')
    else:
        fIn = open(ref)
    with fIn as handle:
        for seq in SeqIO.parse(handle,'fasta'):
            #a list of fragments
            frags = restriction_batch.search(seq.seq)
            locs = []
            for enz_locs in frags.values():
                locs.extend(enz_locs)
            locs = sorted(locs)
            entry = ' '.join(list(map(str, [seq.name] + locs + [len(seq)])))
            yield entry

def digest_handler(reference_fasta,filename_out,enzymes):
    fOut = open(filename_out,'w')
    for entry in digest(enzymes, reference_fasta):
        fOut.write(entry + '\n')


####
####herefore is the old implementation using roll-your-own restriction digests
#### that should be deprecated.

    
    

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
