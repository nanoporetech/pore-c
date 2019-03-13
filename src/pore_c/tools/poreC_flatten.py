import sys 
import re
import uuid

from itertools import combinations, combinations_with_replacement
from collections import defaultdict

class Contact:
    def __init__(self,ch,frag,strand,poss,mapq):
        self.ch = ch
        self.fragID = int(frag)
        if type(strand ) == bool:
            self.strand = 16 if strand else 0
        else:
            self.strand = strand
        self.poss = poss
        self.mapq = mapq
    def __str__(self):
        return ' '.join(list(map(str,(self.ch,self.fragID,self.strand,self.poss,self.mapq))))

class Cwalk:
    def __init__(self,readname=None):
        self.name = readname
#        self.fragIDs = set()
        self.contacts = []

    def from_entry(self,entry):
        l = entry.strip().split()
        L = int((len(l[1:])) / 5)
        self.name = l[0]
        mapqs = l[-L:]
        entries = []

        for x in range(1, 4 * L,4):
            entries.append (l[ x: x + 4])
        
        for x in range(L):
            strand,ch,poss,frag = entries[x]
            if strand == "16":
                strand = True
            else:
                strand = False
            self.add(Contact(ch,frag,strand,poss,mapqs[x]))

    def add(self,monomer):
        self.contacts.append(monomer)
#        self.fragIDs.add(monomer.fragID)


    def sort(self):
        self.contacts.sort(key = lambda c: c.fragID)
        
    def length(self):
        return len(self.contacts)

    def get_chrs(self):
        return tuple([x.ch for x in self.contacts])

    #returns a dictionary of lists of cwalk objects keyed on sorted chromosome sets of size=size
    def flatten(self, size):
        outputs = defaultdict(list)
        chrList = self.get_chrs()
        L = len(self.contacts)

        if len(self.contacts) >= size:
            for c, contact_group in enumerate(combinations(self.contacts, r = size)):
                new_walk = Cwalk("{}_{}".format(self.name,c))
                chrs = tuple(sorted([x.ch for x in contact_group]))
                for contact in contact_group:
                    new_walk.add(contact)
                outputs[chrs].append(new_walk)
            return outputs
        else:
            return None

#        c = 0
#        if len(self.contacts) >= size:
#            for contact_group in combinations(range(L),r = size):
#                new_walk = Cwalk("{}_{}".format(self.name,c))
#                chrs = tuple(sorted([chrList[x] for x in contact_group]))
#                for pos in contact_group:
#                    new_walk.add(self.contacts[pos])
#                outputs[chrs].append(new_walk)
#                c += 1
#            return outputs
#        else:
#            return None
###
    def __str__(self):
        mapString = []
        mapqs = []
        for x in self.contacts:
            tempString = "{strand} {ch} {pos} {frag}".format(strand = x.strand, ch = x.ch, pos = x.poss, frag=x.fragID)
            mapString.append(tempString)
            mapqs.append(x.mapq)

        mapString = ' '.join(mapString)
        quals = ' '.join(list(map(str,mapqs)))
        return "{name} {mappings} {quals}".format(name = self.name, mappings = mapString, quals = quals)



#creates a smart sorted order of chromosomes in the output file
#  handling things like 'chr1,chr2,chr11,chr19,chr20, chr21' and listing chrs with
#  letter names "X,U,W,Y,M..." at the end in alphabetical order
#As written, it sorts in place.
def natural_sort( l ):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

#creates output filehandle
#reads input file
#passes off operations to the functions above as necessary
def flatten_multiway(file_in, file_out, size, sort = True ):
    #gather data by chromosome combination, 
    #   so that it can be printed out in sorted order 
    #   in order to play nicely with juicer

    #this value stores a list of chromosomes visited 
    #   so that they can be written out in an order
    #   that is most compatible with other tools in the
    #   community.
    chr_used = set()
    contact_bundle = defaultdict(list)
    f_out = open(file_out,'w')


    with open(file_in) as f:
        for entry in f:
            walk = Cwalk()
            walk.from_entry(entry)
            if walk.length() < size:
                continue
            if sort:
                walk.sort()
            for contact in walk.contacts:
                chr_used.add(contact.ch)
            flattened_contacts = walk.flatten(size)
            #this is better than a = {**a, **b} because |entries| <<< |chrs| 
            for c , entries in flattened_contacts.items():
                #contacts must be sorted into blocks (e.g. all chr1 intra-chr, then chr1-chr2, chr1-chr3,etc.) which is different from a simple sort.
                for contact in entries:
                    f_out.write(str(contact) + "\n")

    f_out.close()



