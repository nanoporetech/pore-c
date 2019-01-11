import re
import random
from collections import defaultdict
from itertools import combinations, combinations_with_replacement

from pore_c.tools.map_to_frags import Cwalk, Contact

def set_sorted_tuple(data):
    return tuple(sorted(set(data)))

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
def flatten_multiway(file_in, file_out, size, sort ):
    #gather data by chromosome combination, 
    #   so that it can be printed out in sorted order 
    #   in order to play nicely with juicer



    #this value stores a list of chromosomes visited 
    #   so that they can be written out in an order
    #   that is most compatible with other tools in the
    #   community.
    chr_used = set()
    contact_bundle = defaultdict(list)
    with open(file_in) as f:
        for entry in f:
            walk = Cwalk()
            walk.from_entry(entry)
            walk.sort()
            for contact in walk.contacts:
                chr_used.add(contact.ch)
            flattened_contacts = walk.flatten(size)
            #this is better than a = {**a, **b} because |entries| <<< |chrs| 
            for c , entries in flattened_contacts.items():
                contact_bundle[c].extend(entries)

    f_out = open(file_out,'w')

    chr_used = list(chr_used)
    natural_sort(chr_used)

    for chrs in combinations_with_replacement(chr_used,r = size):
        if chrs in contact_bundle:
            for entry in contact_bundle[chrs]:
                f_out.write(str(entry))



