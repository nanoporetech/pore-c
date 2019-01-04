from typing import Pattern, List, NamedTuple, Iterator
from pysam import AlignmentFile
import bisect
import sys

class Contact:
    def __init__(self,ch,frag,strand,poss,mapq):
        self.ch = ch
        self.fragID = frag
        self.strand = 16 if strand else 0
        self.poss = poss
        self.mapq = mapq
    def __str__(self):
        return ' '.join(list(map(str,(self.ch,self.fragID,self.strand,self.poss,self.mapq))))

class Cwalk:
    def __init__(self,readname):
        self.name = readname
        self.fragIDs = set()
        self.contacts = []
    def add(self,monomer):
        if len(self.contacts) == 0 or monomer.fragID != self.contacts[-1].fragID:
            self.contacts.append(monomer)
            self.fragIDs.add(monomer.fragID)
        else:
            #this would be a good point to join fragmented mappings where two alignments fully cover 
            #   a single fragment but were split by the aligner.
            pass

    def sort(self):
        self.sortedContacts = sorted(self.contacts,key = lambda c: c.fragID)

    def __str__(self):
        mapString = []
        mapqs = []
        for x in self.contacts:
            tempString = "{strand} {ch} {pos} {frag}".format(strand = x.strand, ch = x.ch, pos = x.poss, frag=x.fragID)
            if "-1" in tempString:
                continue
            else:
                mapString.append(tempString)
                mapqs.append(x.mapq)

        mapString = ' '.join(mapString)
        quals = ' '.join(mapqs)
        return "{name} {mappings} {quals}".format(name = self.name, mappings = mapString, quals = quals)



#returns ref_frags and ref_IDs
def load_reference(ref_file: str) -> NamedTuple:
    pos = 0
    ref_frags = {}
    ref_IDs = {}
    for entry in open(ref_file):
        l = entry.strip().split()
        l[0] = str(l[0])
        ref_frags[l[0]] = list(map(int,l[1:]))
        ref_IDs[l[0]] = list(range(pos, pos + len(l[1:])))
        pos += len(l[1:])
    return (ref_frags, ref_IDs)


def assign_to_fragment(ref_frags: tuple, ref_IDs: dict, loc: tuple, method: str) -> NamedTuple:

    ch, start, stop = loc

    try:
        assert start < stop
        assert ch in refCoords
        assert ch in refIDs
        assert stop < refCoords[ch][-1] 
    except:
        print('warning: a broken mapping has been seen:\n{}-{} on {} which is {} bp long.'.format(start,stop, ch, refCoords[ch][-1]),file=sys.stderr)

    if method == 'start':
        frag = bisect.bisect_left(refCoords[ch],start + 1 ) -1
        return (refIDs[ch][frag], refCoords[ch][frag])

    if method == 'midpoint':
        frag = bisect.bisect_left(refCoords[ch],int((start+stop)/2) + 1 ) -1
        return (refIDs[ch][frag], refCoords[ch][frag])

    elif method == 'overlap_pct':
        pass


def porec_iterator(input_bam: str):
    aligns = []
    current_read_name = None
    for align in bam:
        if current_read_name is None:
            current_read_name = align.query_name
            aligns.append(align)
        elif current_read_name == align.query_name:
            aligns.append(align)
        else:
            yield (current_read_name, sorted(aligns,key=lambda x: x.query_alignment_start))
            aligns = [align]
            current_read_name = align.query_name
    yield (current_read_name, sorted(aligns,key=lambda x: x.query_alignment_start))


def map_to_fragments(input_bam: str, ref_file: str, output_file: str, method: str) -> None:
    
    ref_frags, ref_IDs = load_reference(ref_file)

    f_out = open(output_file, 'w')

    for name, aligns in porec_iterator(input_bam):
        walk = Cwalk(name)
        for monomer in aligns:
            loc = (monomer.reference_name, monomer.reference_start,monomer.reference_end)
            frag,pos = assign_to_fragment(ref_frags, ref_IDs, loc, method = method)
            new_contact = Contact(monomer.reference_name, frag, monomer.is_reverse, pos, monomer.mapping_quality)
            walk.add(new_contact)
        f_out.write(str(walk))
    
    f_out.close()

            
