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
reSeq = re.compile(str(sys.argv[2]).encode('utf-8'))

for entry in SeqIO.parse(sys.argv[1],"fasta"):
    refOut = [entry.name]
    print('reading {}, of {} bp.'.format(entry.name,len(entry.seq)),file=sys.stderr)
    #case where sequence contains no sites
    seq = str(entry.seq).encode('utf-8')
    if not re.search(reSeq,seq):
        refOut.append(str(len(seq)))
        print(' '.join(refOut))
        continue
    #iterate over sites
    for site in re.finditer(reSeq,seq):
        refOut.append(str(site.start()))
    print(' '.join(refOut))
    
