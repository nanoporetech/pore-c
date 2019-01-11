#!/usr/bin/env python

import random
import re
import argparse
import sys
from itertools import combinations,combinations_with_replacement
from collections import defaultdict

#takes in a list, creates a set, sorts the entries and makes it hashable (as a tuple)
def setSortedTuple(data):
    return tuple(sorted(set(data)))

def chrSmartsort(chrs):
    strChrs = set()
    intChrs = set()
    complexIntChrs = set()
    for entry in chrs:
        #scan for simple int chrs
        try:
            assert int(entry)
            intChrs.add(int(entry))
        except:
            strChrs.add(entry)
    l = list(strChrs)
    for entry in l:
        # scan for complex int chrs, e.g., chr1, chr2, etc.
        prefix_re = re.search('^(?P<prefix>\D+)[0-9]+',entry)
        if prefix_re:
            prefix = prefix_re.group('prefix')
            try:
                newEntry = int(entry.replace(prefix,''))
                if newEntry not in intChrs:
                    intChrs.add(newEntry)
                    complexIntChrs.add(newEntry)
                    strChrs.remove(entry)
                else:
                    raise Error('Duplicate detected {} and {}'.format(entry, newEntry))
            except:
                pass
    output = []
    for entry in sorted(intChrs):
        if entry in complexIntChrs:
            output.append(prefix+str(entry))
        else:
            output.append(str(entry))
    for entry in sorted(strChrs):
        output.append(entry)
    return output
                        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "This tool takes in a poreC and decomposes the multiway contacts to a specified dimension. Multiway contacts of length less than the specified dimensions are ignored.")
    parser.add_argument('--poreC' ,'-p', help="A pore-C file to be decomposed into constituent contacts.")
    parser.add_argument('--sort' ,'-s', action="store_true", help="When printing the contacts for a single read, sort the contacts by the fragmentID")
    parser.add_argument('--setSize' , type = int, default = 2, help="The size of fragment sets outputted. by default, the output is pairwise data.")

    args = parser.parse_args()

    #this mechanism is included so that the contact map is printed out in an orderly fashion. first chr1 intracontacts, then chr1,chr2  chr1,chr3... chr2,chr3 ... chr3,chr4...
    chrsUsed = set()

    chrData = defaultdict(list)
    counts = 0
    print("gathering chromosome data from {}".format(args.poreC), file = sys.stderr)
    for line in open(args.poreC):
        l = line.strip().split()
        readID = l[0]
        L = len(l[1:])
        readCount = int(L / 5)
        #lines with fewer than the specified size of the desired fragment
        #  are discarded.        
        if readCount < args.setSize:
            continue

        entries = []
        #collect the data by pair of chromosomes
        for pos in range(1,readCount*4,4):
            entry = l[pos:pos+4]
            entries.append("\t".join(entry))
            if entry[1] not in chrsUsed:
                chrsUsed.add(entry[1])
        mapqs = l[-readCount:]
        assert len(mapqs)== len(entries)
        mapqDict = {}
        for pos in range(len(mapqs)):
            mapqDict[entries[pos]] = mapqs[pos]
        for _comb in combinations(entries,r=args.setSize):
            if args.sort:
                comb = sorted(_comb, key=lambda x: int(x.split()[3]))
            else:
                comb = _comb
            #these should already be implicitly sorted as per fragmentIDs
            # IF sorted=True for output...THIS IS A BUG.
            chrs = tuple(map(lambda x: x.split()[1],comb))
            chrData[chrs].append((comb,mapqDict))
            counts += 1

    chrsUsed = chrSmartsort(list(chrsUsed))

    a = list(map(lambda x: (x[0], len(x[1])), chrData.items()))
    print("exhibit A:",a,file=sys.stderr)
                    
    counts = 0
    for chrs in combinations_with_replacement(chrsUsed,r = args.setSize):
        print('printing out {} combinations from chromosomes {}'.format(len(chrData[chrs]),' '.join(chrs)), file = sys.stderr)
        for comb , mapqDict in chrData[chrs]:
            #each contact is given a unique contact identifier composed of the poreC readID and a unique number.
            #   all decomposed contacts originating from the same poreC read should therefore be identifiable by the
            #   readID component of their unique contact identifier.

            print('{readID}_{count}\t{contacts}\t{mapqs}'.format(readID=readID, count = counts, contacts = '\t'.join(comb),mapqs = '\t'.join(map(lambda x: str(mapqDict[x]), comb))))        
            counts += 1

    print('{} {}-way contacts printed.'.format(counts,args.setSize),file=sys.stderr)
