#!/usr/bin/env python

#turns a .hic + .hicREF file to a .bed based sparse representation 
#  used by a variety of non-erez programs for analysis

import sys
from collections import Counter

ref = sys.argv[2]

fragToChr = {}
fragToInterval = {}
curFrag = 0
for entry in open(ref):
    l = entry.strip().split()
    fragToInterval[curFrag] = (0,int(l[1]))
    curFrag += 1
    for pos in range(2, len(l) -1 ):
        fragToChr[curFrag] = l[0]
        fragToInterval[curFrag] = (l[pos], l[pos-1])
        #prints out the reference header information:
        #chr start end binID
        print('{}\t{}\t{}\t{}'.format(l[0],l[pos],l[pos-1]+1,curFrag))
        curFrag += 1

contacts = Counter()

for entry in open(sys.argv[1]):
     l = entry.strip().split()
     coord = tuple(sorted( [int(l[4]) , int(l[8])]))
     contacts[coord] += 1


for coord, count in contacts.items():
    #prints out the data as:
    #bin1 bin2 count
    print("{}\t{}\t{}".format(coord[0],coord[1], count))
