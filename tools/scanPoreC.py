#!/usr/bin/env python

### Takes in:
###   a .fragment file (one fragmentID per line) and 
###   a .poreC file

### returns:
###  any poreC entry containing any of the fragments on
###  the list.

import sys

fragsIn = sys.argv[1]

poreCin = sys.argv[2]

frags = set()

for entry in open(fragsIn):
    frags.add(entry.strip())

for entry in open(poreCin):
    printFlag = False
    l = entry.strip().split()[1:]
    z = len(l)
    if z %  5 != 0:
       raise ValueError('oops. bad poreC entrty.')

    fragNum = int(z / 5)
    #print(entry,l,z, fragNum)
    mapqs = l[-fragNum:]
    entries = []
    for pos in range(0, 4 * fragNum, 4):
        #print ("entry:",tuple(l[pos:pos+4]))
        entries.append(tuple(l[pos:pos+4]))
    #print("entries:",entries)
    for monomer in entries:
        if monomer[3] in frags:
            printFlag = True
    if printFlag:
        print(entry.strip())
        
