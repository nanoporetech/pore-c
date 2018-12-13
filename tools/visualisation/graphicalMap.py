#!/usr/bin/env python 

import sys
import hashlib
import pysam
import re
import matplotlib
matplotlib.use('agg')
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle 

locus = re.compile("(?P<ch>([0-9]+|X|Y|M)):(?P<start>[0-9]+)-(?P<end>[0-9]+):(?P<strand>(-1|1))")

labels = ".".join(sys.argv[1].split('.')[:-1])
print(labels)

fIn = pysam.AlignmentFile(sys.argv[1])
#includedAlns = pysam.AlignmentFile(sys.argv[2])
if len(sys.argv) > 2:
    groundTruth = open(sys.argv[2])

    GT = defaultdict(list)
    sizes = {}
    #creates an ordered list for each read indicating the length of the fragments
    #  this will be drawn onto the graphs to indicate where fragments actually
    #  should be ending/beginning
    for line in groundTruth:
        l = line.strip().split()
        name = "_".join(map(str, map(int, l[0].split("/"))))
        #what is sizes? The number of 
        sizes[name] = int(name.split('_')[-1])
        for x in re.finditer(locus,line):
            if name not in GT:
                GT[name].append(0)
            GT[name].append(int(x.group("end")) - int(x.group("start")) + GT[name][-1])

else:
    GT = defaultdict(list)
    sizes = 0

#alns = set()
#for x in includedAlns:
#    alns.add(hashlib.sha1(str(x).encode('utf-8')).hexdigest())

mult = 100
lastRead = False
for x in fIn:
    currentRead = x.query_name#"_".join(map(str, map(int, x.qname.split("/"))))
    if lastRead == False: #case: first read, setting things up
        fig = plt.figure()
        ax = fig.subplots(1)
        layer = 0
        maxlen = 0
        lastRead = currentRead

    if lastRead != currentRead: #a read set has closed
        #save plot, then clear the slate
        fn = "{}.{}.png".format(labels,lastRead).replace("/","_")
        print("plotting {}. {} alignments found".format(fn, int( layer / mult)), file=sys.stderr)
        ax.get_yaxis().set_visible("off")
#        if sizes != 0:
#            ax.set_xlim(left = 0, right = sizes[lastRead] )
#        else:
        ax.set_xlim(left = 0 , right = pos * 1.05)
        ax.set_ylim(bottom = -20, top = layer + mult)

        #ground truth dashed lines
        if currentRead in GT:
            ax.vlines(GT[currentRead],0,layer + mult, linestyles="dotted")
        [line.set_marker('None') for line in ax.get_yticklines()] 
        ax.get_yaxis().set_visible(False)
        ax.set_yticklabels([])
        print(fn)
        fig.savefig(fn,dpi=1200)
        plt.close(fig)

        #new plot initialised
        fig = plt.figure()
        ax = fig.subplots(1)
        lastRead = currentRead
        layer = 0
        maxlen = 0

    #plot with rectangles plotted along a line
    #different rectangles correspond to different cigar operations
    if x.is_reverse:
#        cigar = x.cigar
        cigar = x.cigar[::-1]
    else:
        #sometimes a cigar is just a cigar
        cigar = x.cigar
    pos = 0
    aligned_length = 0
#    print(cigar)
    for op, length in cigar:
        #        print(op,length,pos,layer)
        if op == 0: #match
            rect = Rectangle((pos,layer),length, mult / 2, facecolor = "chartreuse")
        elif op == 2: #insertion (or deletion in reference)
#            ax.plot([pos,pos ], [layer,layer+mult / 2], color="violet", linewidth = 1)
            rect = Rectangle((pos,layer-10), 1, mult / 5, facecolor = "violet")
#            continue
        elif op == 1: #deletion (or insertion in reference)
            rect = Rectangle((pos,layer),length, mult / 5, facecolor = "darkred")
        elif op == 3: #skipped region from reference
            rect = Rectangle((pos,layer),length, mult / 2.2, facecolor = "grey")
        elif op == 4: #soft clip
            rect = Rectangle((pos,layer),length, mult / 3, facecolor = "lavender")
        elif op == 5: #hard clip
            rect = Rectangle((pos,layer),length, mult / 3 , facecolor = "silver")
        elif op == 6: #padding
            rect = Rectangle((pos,layer),length, mult / 3, facecolor = "grey")
        elif op == 7: #match
            rect = Rectangle((pos,layer),length, mult / 2, facecolor = "chartreuse")
        elif op == 8: #mismatch
            rect = Rectangle((pos,layer),length, mult / 2.1, facecolor = "darkgreen")
        if op != 2:
            pos += length
        ax.add_patch(rect)
   
    #if pos > sizes[currentRead]:
    #    print('warning. Read {} of cigar length {} exceeded read length of {}.'.format(currentRead,sizes[currentRead],pos))
    #    print(x)
    if x.is_reverse:
        ax.plot(pos - 10, layer - 10,"k<")
    else:
        ax.plot(10, layer - 10 ,"k>")

        #readstartrect = Rectangle((0,layer + layer / 3),100, mult / 3, facecolor = "black")

    layer += mult 
