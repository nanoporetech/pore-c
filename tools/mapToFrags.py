
#!/usr/bin/env python

## Because the notion of order becomes important to concatemers, handling of order and order printing
#  has to be done with some flexibility at first in order to accomodate several potential downstream
#  applications, e.g., some operations may be easier with sorted data, but original order may also be
#  important, so this script will be written so that it can output either one at the user's discretion.

class contact:
    def __init__(self,ch,frag,strand,poss,mapq):
        self.ch = ch
        self.fragID = frag
        #strand is from pysam.Alignment.is_reverse, using the pre .hic specifications
        self.strand = 16 if strand else 0
        self.poss = poss
        self.mapq = mapq
    def __str__(self):
        return ' '.join(list(map(str,(self.ch,self.fragID,self.strand,self.poss,self.mapq))))

    # all comparisons eval by fragmentID
    def __lt__(self, other):
        if self.fragID < other.fragID:
            return True
        return False
    
    def __gt__(self, other):
        if self.fragID > other.fragID:
            return True
        return False

    def __le__(self, other):
        if self.fragID <= other.fragID:
            return True
        return False
    
    def __ge__(self, other):
        if self.fragID >= other.fragID:
            return True
        return False        


class cwalk:
    def __init__(self,readname):

        self.name = readname
        self.fragIDs = set()
        self.contacts = []

    def add(self,monomer):
        #this addresses cases where the alignment of a fragment is broken up
        #(i.e., two sequential mappings are to the same fragment)
        if len(self.contacts) == 0 or  monomer.fragID != self.contacts[-1].fragID:
            self.contacts.append(monomer)
            self.fragIDs.add(monomer.fragID)
        
        
    def boop(self):
        self.contacts.sort(key=lambda c: c.fragID)
        
    def decompose(self, setSize, sorted=False):
        ###Nb. This has not been tested and is superceded by the pore-c porecflatten.py script.
        self.length = len(self.contacts)
        poss = list(range(self.length))
        results = []
        uniqueID = 0
        for entries in combinations(self.contacts(),r=setSize):
            mappings = []
            quals = []
            # create a sorted set of entries
            sortedEntries = sorted(entries, lambda x: x.fragID)

            #this may be unnecessary if combinations works correctly
            tse = tuple(sortedEntries)
            if tse in sortedEntrySet:
                continue
            sortedEntrySet.add(tse)

            if sorted:
                Z = sortedEntries
            else:
                Z = entries
            
            for entry in Z:
                tempString = "{strand} {ch} {pos} {frag}".format(strand=entry.strand,ch = entry.ch, pos = entry.poss,frag = entry.frag)
                if "-1" in tempString:
                    continue
                else:
                    mappings.append(tempString)

            quals = ' '.join(quals)
            mappings = ' '.join(mappings)
            results.append('{readID}-{UID} {mappings} {quals}\n'.format(readID = self.name, UID = uniqueID, mappings = mappings, quals = quals))
            uniqueID += 1
        return ''.join(results)
    
    #by default, this returns the contacts in the order seen in the molecule
    # if sorted == True, then it sorts the contacts before returning them
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

def assignToFragment(refCoords,refIDs,loc,method="start",pct = 0.7):

    ch,start,stop = loc

    assert start < stop
    assert ch in refCoords
    assert ch in refIDs
    
    try:
        assert stop < refCoords[ch][-1]
    except:
        pass

    #Hi-C data requires a map-filtering step based on the relationship between an alignment
    # and expected fragment boundaries. the "start" method here only assigns a mapping to a 
    # fragment based on the start position of the mapped interval. If a filtering step were
    # to be implemented based on the quality of the mapping of a pore-C monomer to a given
    # fragment (either based on pct aligned, alignment quality, etc.) this would be where it
    # is implemented.    
    if method == 'start':
        frag = bisect.bisect_left(refCoords[ch],start +1) - 1
        return refIDs[ch][frag],refCoords[ch][frag]
    elif method == 'overlap_pct':
        pass
    elif method == 'all':
        raise ValueError("this method isn't supported yet.")
    
        firstFrag = bisect.bisect_left(fragRef[ch],start) 
        lastFrag = firstFrag + bisect.bisect_left(fragRef[ch][firstFrag:],stop) - 1
        return list(range(firstFrag,lastFrag+1))
    
if __name__ == "__main__":

    import argparse
    import pysam as ps
    #TODO [ ]: add paf file support in addition to .bam support. 
    #          just create a file handling obj for paf or bam and iterate over that. make sure paf handler has same attribute names as sam.
    import bisect
    import sys
    
    parser = argparse.ArgumentParser(description = "This tool takes in a bam file and fragment map, maps each alignment from each read to a fragment and returns a poreC file which can be further manipulated and explored with X Y Z...tools.")
    parser.add_argument('--alignments' ,'-a', help="A bam or paf alignment file of pore-C data to a genome of interest.")
    parser.add_argument('--fragmentref' ,'-f', help="The fragment reference file, in the style of the juicer reference file format.")
    parser.add_argument('--setSize' , default = 1000, help="The size of fragment sets outputted. by default, the output is effectively full length concatemers. This has not ")
    parser.add_argument('--sort' ,'-s', type=bool, default=False, help="When printing the contacts for a single read, sort the contacts by the fragmentID.")
    parser.add_argument('--method' ,'-m',default ="start", help="The assignment method. Accepted values are 'start','overlap_pct','all',and maybe some others. Default is 'start'.")

    args = parser.parse_args()

    #load reference fragment file
    refFrags = {}
    refIDs = {}
    pos = 0
    for entry in open(args.fragmentref):
        l = entry.strip().split()
        l[0] = str(l[0])
        refFrags[l[0]] = list(map(int,l[1:]))
        refIDs[l[0]] = list(range(pos, pos + len(l[1:])))
        pos += len(l[1:])
        print('ref {}: {} ... {}'.format(l[0], ' '.join(list(map(str,refIDs[l[0]][:3]))), ' '.join(list(map(str,refIDs[l[0]][-3:])))),file=sys.stderr)

    print('Loaded reference genome {} with {} total fragments.'.format(args.fragmentref,pos),file=sys.stderr)

    #a read of interest, for testing purposes to be sure all fragments are being properly captured and in correct order
    ##V##ROI = "0000be7c-f480-4cd8-9326-b8e16a3ade26"
    ##V##print ('looking at ROI: {}'.format(ROI),file=sys.stderr)
    #bamFile/pafFile
    lastRead = False
    count = 0
    print("Mapping {} to fragments...".format(args.alignments),file = sys.stderr)
    #####HERE BEGINS THE SPONGE!!!
    for entry in ps.AlignmentFile(args.alignments):
        
        ##V##if entry.query_name == ROI:
        ##V##    print (">>>",entry.reference_name, entry.reference_start,entry.reference_end)

        if entry.is_unmapped:
            if lastRead:
                pass
            else:
                lastRead = entry.query_name
                entries = []
            continue
        
        if not lastRead:
            entries = [entry]        
        elif entry.query_name != lastRead:
            try:
                assert entries
            except:
                print(entry.query_name,lastRead,entry.is_unmapped)
            if len(entries) > 1:
                #reads must be ordered according to location in the poreC read
                entries = sorted(entries, key=lambda x: x.query_alignment_start)
                walk = cwalk(lastRead)
                for monomer in entries:
                    frag,poss = assignToFragment(refFrags,refIDs,(monomer.reference_name,monomer.reference_start,monomer.reference_end),method=args.method)
                    new_contact = contact(monomer.reference_name,frag,monomer.is_reverse,poss,str(monomer.mapq))
                    walk.add(new_contact)
                ##V##if lastRead == ROI:
                ##V##    print(walk)
                ##V##    for monomer in entries:
                ##V##        print("sorted:",monomer.reference_name, monomer.reference_start,monomer.reference_end)

                #print(walk)
                print(walk)

            entries = [entry]
            
            count +=1
            ##V##if lastRead == ROI:
            ##V##    print([x.fragID for x in walk.contacts])
            ##V##    print(walk) #for compatibility elsewhere
            
            #fragments - this doesn't yet support "all" as handling of multi-fragments isn't in place yet.

        elif entry.query_name == lastRead:
            entries.append(entry)

        lastRead = entry.query_name

    if not lastRead:
        for monomer in entries:
            frag,poss = assignToFragment(refFrags,refIDs,(monomer.reference_name,monomer.reference_start,monomer.reference_end),method=args.method)
            new_contact = contact(monomer.reference_name,frag,monomer.is_reverse,poss,str(monomer.mapping_quality))          
            walk.add(new_contact)            
            if args.sort == True:
                walk.boop()
            if len(walk.contacts) > 1:        
                #print(count, walk)
                print(walk) #compatible with other tools elsewhere without the index line
                count += 1
    else:
        print('no mappings identified.',file=sys.stderr)
    


