#!/usr/bin/env python

#takes in a reference genome and generates the necessary reference
#   map used by pre for inputting our data into a pre-.hic format.
"""
format:
chr site1 site2 site3 ... chrLength

as in:
1 11160 12411 12461 ... 249250621
2 11514 11874 12160 ... 243199373
3 60138 60662 60788 ... 198022430
"""
import sys
import re
from Bio import SeqIO

#the restriction enzyme motif sequence, as a regex, if non-palindromic
#reSeq = re.compile(str(sys.argv[2]).encode('utf-8'))

###
rcDict = dict(zip('ACGTWSMKRYBDHVNacgtwsmkrybdhvn-','TGCAWSKMYRVHDBNtgcawskmyrvhdbn-'))

def revcomp(seq):
    """Return the reverse complement of a string:

    :param seq: DNA Sequence to be reverse complemented
    :type seq: str
    :returns: A new string that is the reverse complement of seq
    :rtype: string

    """
    return ''.join([rcDict[nuc] for nuc in  seq[::-1]])

def REGEXsite(seq):
    s = seq.replace("N","[ACGT]")
    s = s.replace("V","[ACG]").replace("H","[ACT]").replace("D","[AGT]").replace("B","[CGT]")
    s = s.replace("W","[AT]").replace("S","[CG]").replace("M","[AC]").replace("K","[GT]").replace("R","[AG]").replace("Y","[CT]")
    return s

site_raw = sys.argv[2]
print('raw site:',site_raw,file = sys.stderr)
fwd_rev_site = "({}|{})".format(site_raw,revcomp(site_raw))
print('fwdrev site:',fwd_rev_site,file = sys.stderr)
re_site = REGEXsite(fwd_rev_site)
print('REGEX site:',re_site,file = sys.stderr)
site_compiled = re.compile(re_site)

fragMaps = {}
#keyed on a (ind,len) tuple so that chrs can be sorted by length (getting around non-alphanumeric sorting problems)
# while still guaranteeing key uniqueness by the indicator

###

#accepts gzip compressed reference fasta
if '.gz' in sys.argv[1]:
    fIn = gzip.open(sys.argv[1], 'rt')
else:
    fIn = open(sys.argv[1])

ind = 0
with fIn as handle:
    for entry in SeqIO.parse(sys.argv[1],"fasta"):
        refOut = [entry.name]
        print('reading {}, of {} bp.'.format(entry.name,len(entry.seq)),file=sys.stderr)
        #case where sequence contains no sites
        seq = str(entry.seq)
        if not re.search(site_compiled,seq):
            print('no {} sites found in {}.'.format(re_site,entry.name),file=sys.stderr)
            refOut.append(str(len(seq)))
            print(' '.join(refOut))
            continue
        #iterate over sites
        for site in re.finditer(site_compiled,seq):
            refOut.append(str(site.start()))
        #add the seq length to cap off the entry
        refOut.append(str(len(seq)))
        fragMaps[(ind,len(seq))] = ' '.join(refOut)
    #    print(' '.join(refOut))
        ind += 1


fragMaps = sorted(list(fragMaps.items()),key = lambda x: x[0][1], reverse = True)

for key,value in fragMaps:
    print(value)
