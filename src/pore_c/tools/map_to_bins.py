from typing import Pattern, List, NamedTuple, Iterator, Union, Dict, Tuple

import gzip
from pybedtools import BedTool
from dataclasses import dataclass
from collections import Counter

@dataclass
class HicContact(object):
    read_id: str
    strand1: bool
    chr1: str
    pos1: int
    frag1: int
    strand2: bool
    chr2: str
    pos2: int
    frag2: int
    mapq1: int
    mapq2: int

    @classmethod
    def from_hictxt(cls, hictxt):
        vals = hictxt.strip().split(' ')
        assert len(vals) == 11
        for field in [3,4,7,8,9,10]:
            vals[field] = int(vals[field])
        if vals[1] == "16":
            vals[1] = True 
        else:
            vals[1] = False
        if vals[5] == "16":
            vals[5] = True 
        else:
            vals[5] = False
        return cls(*vals)

    def to_midpoint_bed_pair(self, fragment_reference):

        ref_ch1, midpoint1 = fragment_reference[self.frag1]
        ref_ch2, midpoint2 = fragment_reference[self.frag2]
        assert ref_ch1 == self.chr1
        assert ref_ch2 == self.chr2
        
        data_template = "{ch}\t{st}\t{en}\t{read_id}\n"

        return data_template.format(ch = ref_ch1,st=midpoint1, en = midpoint1+1, read_id = self.read_id),\
            data_template.format(ch = ref_ch2,st=midpoint2, en = midpoint2+1, read_id = self.read_id)


def bin_hic_data(input_hictxt: str, output_bin_matrix: str, fragment_reference: str, bin_reference: str) -> None:
    frag_midpoint_reference = {}
    with gzip.open(fragment_reference) as handle:
        for entry in handle:
            l = entry.decode('utf-8').strip().split()
            frag_midpoint_reference[int(l[3])] = (l[0], (int(l[1]) + int(l[2]) // 2))

    #gather contacts into a list of strings to be concatenated 
    # . and converted into a single bedtool object for mashing against the reference bin bedfile
    contacts = []
    for entry in map( HicContact.from_hictxt, open(input_hictxt)):
        c1, c2 = entry.to_midpoint_bed_pair(frag_midpoint_reference)
        contacts.append(c1)
        contacts.append(c2)

    contact_bed = BedTool("\n".join(contacts),from_string= True)
    
    bins_bed = BedTool(bin_reference)
    map_bed = contact_bed.intersect(bins_bed,wao=True) #should this be outputted as a helper file

    contacts_seen = {}
    sparse_matrix = Counter()
    for entry in map_bed:
        l = str(entry).strip().split()
        if l[3] in contacts_seen:
            #I think sort order is preserved here by virtue of the pipeline, though it may be worthwhile to
            # . explicitly enforce it
            if contacts_seen[l[3]] > int(l[7]):
                x1, x2 = int(l[7]) , contacts_seen[l[3]]
            else:
                x1, x2 = contacts_seen[l[3]], int(l[7])
            sparse_matrix[(x1,x2)] += 1
        else:
            contacts_seen[l[3]] = int(l[7])

    print(sparse_matrix)

    f_out = open(output_bin_matrix,'w')
                        

    for coords, val in sparse_matrix.items():
        x,y = coords
        f_out.write("{} {} {}\n".format(x,y,val))

    f_out.close()

