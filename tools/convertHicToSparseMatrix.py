
import struct

import gzip
import straw
import sys
from dataclasses import dataclass

@dataclass
class bedEntry:
    ch: str
    start: int
    end: int
    frag_id: int

    @classmethod
    def from_string(cls,entry):
        l = entry.strip().split('\t')
        return cls( *l)

###########
#an excerpt from the read_hic_header.py file solely for extracting the available chromosome names
###########
def readcstr(f):
    buf = ""
    while True:
        b = f.read(1)
        b = b.decode('utf-8', 'backslashreplace')
        if b is None or b == '\0':
            return str(buf)
        else:
            buf = buf + b

def get_chrs(req):
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        print('This does not appear to be a HiC file magic string is incorrect')
        sys.exit(1)
    version = struct.unpack('<i',req.read(4))[0]
#    print('HiC version:')
#    print('  {0}'.format(str(version))) 
    masterindex = struct.unpack('<q',req.read(8))[0]
#    print('Master index:')
#    print('  {0}'.format(str(masterindex)))
    genome = ""
    c=req.read(1).decode("utf-8") 
    while (c != '\0'):
        genome += c
        c=req.read(1).decode("utf-8") 
#    print('Genome ID:')
#    print('  {0}'.format(str(genome))) 
    # read and throw away attribute dictionary (stats+graphs)
#    if (verbose is 1):
#        print('Attribute dictionary:')
    nattributes = struct.unpack('<i',req.read(4))[0]
    for x in range(0, nattributes):
        key = readcstr(req)
        value = readcstr(req)
#    if (verbose is 1):
#        print('   Key:{0}'.format(key))
#        print('   Value:{0}'.format(value))
    nChrs = struct.unpack('<i',req.read(4))[0]
#    print("Chromosomes: ")
    names = []
    lengths = []
    for x in range(0, nChrs):
        names.append(readcstr(req))
        lengths.append(struct.unpack('<i',req.read(4))[0])

    nBpRes = struct.unpack('<i',req.read(4))[0]
    #print("Base pair-delimited resolutions: ")
    bin_resolutions = []
    for x in range(0, nBpRes):
        bin_resolutions.append(struct.unpack('<i',req.read(4))[0])
        #print('   {0}'.format(res))
    nFrag = struct.unpack('<i',req.read(4))[0]
    fragment_resolutions = []
    #print("Fragment-delimited resolutions: ")
    for x in range(0, nFrag):
        fragment_resolutions.append( struct.unpack('<i',req.read(4))[0])
    return names,lengths,bin_resolutions,fragment_resolutions


hic_file_in = open(sys.argv[1],'rb')
resolution = int(sys.argv[2])


names,lengths,bin_resolutions,fragment_resolutions = get_chrs(hic_file_in)

if resolution not in bin_resolutions:
    raise ValueError("Resolution {} is not supported by this .hic file. Choose from among the following bin level resolutions:\n{}\nor the following fragement level resolutions:\n{}\n".format(resolution,' '.join(list(map(str,bin_resolutions))),' '.join(list(map(str,fragment_resolutions)))))

reference_bed = sys.argv[3] #this needs to be the digest reference bed file

frags = {}
chrs = set()
for entry in map(bedEntry.from_string, gzip.open(reference_bed,'rt')):
    if "GL" not in entry.ch:
        chrs.add(entry.ch)
    frags[(entry.ch, int(entry.start))] = entry.frag_id


#print(list(frags.items()))
#disregard the low resolution matrix view from the dump
names.remove("All")

print("hicfile:",set(names))
print("bedfile:",set(chrs))

if not set(names).issubset( chrs):
    print( set(names) - chrs)
    print(chrs - set(names))
    raise ValueError("The reference file does not have all the chromosomes described in the .hic binary file.")

matrix_out = open( sys.argv[4],'w')

exp_lens = dict(zip(lengths,names))

for x in range(len(names)):
    for y in range(x,len(names)):
        ch1 = names[x]
        ch2 = names[y]
        if ch1 == "All" or ch2 == "All":
            continue
        print(ch1,ch2)
        bin_loc1,bin_loc2,counts = straw.straw("NONE",sys.argv[1],ch1,ch2,"BP", resolution)
        for pos in range(len(counts)):
            bin1 = frags[(ch1,bin_loc1[pos])]
            bin2 = frags[(ch2,bin_loc2[pos])]


            matrix_out.write("{bin1}\t{bin2}\t{count}\n".format(bin1= bin1,bin2=bin2,count = counts[pos]))


