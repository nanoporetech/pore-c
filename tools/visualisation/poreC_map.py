#!/usr/bin/env python


#takes in a poreC file, a resolution, and a reference genome 

#generates a plot of all indicated poreC reads along the reference

#plotting specifications will be such that poreC monomers will be 
# will be plotted as rectangles linked by a line with the color
# indicating the order of the mappings. Rectangles will either be below
# or above a poreC guideline, below if they're minus strand sequences
# and above if they're plus strand. Because the resolution of a rectangle
# of fragment size in a genome is not possible (or even of a bin sized
# demarcation), no translation from fragment to bin to mapping location
# will be done. Instead, poreC monomers will be mapped directly to their location (as opposed to
# fragment or bin).

# For assistance, bin demarcations might be added in the future for higher resolution plots.

# Note: Because the reference is dictated by the refFile entries, it's possible
# to plot just the mappings to a specific chromosome or subset of chromosomes
# by generating a reference with only the chromosomes of interest in the order
# of interest. In that case, a flag will be added which can enable the addition of
# guidelines at a specified resolution (1M or 10M guidelines, for example).

import sys

import numpy as np

from collections import defaultdict

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.patches import Rectangle
poreC_file = sys.argv[1]
prefix = poreC_file.replace('.poreC','')

refFile = sys.argv[2]

#this will only be used in order to add minor guidelines for bin boundaries
#binSize = int(sys.argv[3])

###first, load the reference
###   in this case, this simply means grabbing the reference lengths
###   and fragment lengths by number (for plotting purposes)
print('loading reference genome data.',file=sys.stderr)

curPos = 0
chrList = []
chrPosList_ordered = []
curFrag = 0
fragSize = {} #stores an ordered list of fragment sizes for plotting rectangles
chrPosList = {} #stores the cumulative start position of each Chr in the genome
for entry in open(refFile):

    l = entry.strip().split()
    print("loading {}".format(l[0]), file=sys.stderr)
    chrList.append(l[0])
    chrPosList[l[0]] = curPos
    chrPosList_ordered.append(curPos)
    fragSize[curFrag] = int(l[1])
    fragID = np.arange(curFrag,curFrag + len(l[1:]))
    starts = np.array([0] + list(map(int,l[1:-1])))
    stops = np.array(list(map(int,l[1:])))
    diffs = stops - starts
    for pos in range(1, len(starts)):
        fragSize[fragID[pos]] = diffs[pos]
    curPos += int(l[-1])
    curFrag += len(l[1:])

###load the poreC file, storing a sequence of locations for each fragment set
#each item in poreCData will be a list of read (ch, locations). the ultimate loc
# of each mapping's point will be (strand, chrPosList[ch] + locations)
poreCData = []
for entry in open(poreC_file):
    l = entry.strip().split()
    if len(l[1:]) % 5 != 0:
        raise ValueError('bad entry.')
    mapCount = int(len(l[1:]) / 5)
    
    mappings = []
    for pos in range(1,mapCount*4, 4):
        mappings.append((int(l[pos]), chrPosList[l[pos+1]] + int(l[pos+2])))
    poreCData.append(mappings)

plt.clf()
fig,ax = plt.subplots()
spacing = 10
layer = 1
readCount = 0
#cw = cm.get_cmap('jet')
for readNum, reads in enumerate(poreCData):
    if len(reads) < 2:
        continue
    readCount += 1
    minRead = min(list(map(lambda x: x[1], reads)))
    maxRead =  max(list(map(lambda x: x[1], reads)))
    #plot the guideline from the leftmost location to the rightmost location, dashed and light in ONT color scheme
    ax.hlines(layer,minRead,maxRead, linestyle = "-.", linewidth = .1, color='#357BA1')

    print([minRead, maxRead], [layer,layer])
#    plt.plot([minRead, maxRead], [layer,layer] ,linestyle = "-.", linewidth = 10 , color = 'k')#'#357BA1')
    L = len(reads)
    #set color for all reads
    colorvals = [float(x / (L-1)) for x in range(int(L))]
    print(colorvals)
    for col, mapping in enumerate(reads):
        if mapping[0] == 16:

            offset = -spacing 
        else:
            #the vertical size should be spacing / 2 regardless of plus or minus strand
            offset = 0
        print('+++read+++')
        print('strand:', mapping[0])
        print("x,y:", (mapping[1], layer))
        print("width:", fragSize[readNum])
        print("height:",spacing)
        print("colour:",colorvals[col])
        rect = Rectangle((mapping[1],layer),fragSize[readNum],spacing/2,color=cm.jet(colorvals[col]))
        ax.add_patch(rect)
    layer += spacing
        
###plot

#plt.imshow(contactMatrix, norm=colors.LogNorm(vmin=1, vmax=contactMatrix.max()), cmap=plt.get_cmap("gray_r"))
ax.vlines(chrPosList_ordered,0,readCount * spacing, linestyle = ":", linewidth = .1, alpha=0.3, color='r')# color = '#357BA1')

#ax.set_xticks()
#ax.set_xticklabels()
ax.set_xticks(chrPosList_ordered)
ax.set_xticklabels(chrList)
ax.set_ylim(0,readCount * spacing)

#add colorbar
#cbar = plt.colorbar(orientation="horizonal",shrink=.1)
plt.tick_params(axis="x", which='major', labelleft=False, left=False, labelsize=4)
plt.xticks()
#########TODO ADD A COLORSCHEME LEGEND
#   READ------------->BLUE
#           READ
plt.savefig('{}.poreC_mappings.png'.format(prefix),dpi=2000)



