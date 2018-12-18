#!/usr/bin/env python


#takes in a hic file, a resolution, and a reference genome 

#generates a matrix in the style of the poreC.py script

# simple as that

import sys

import numpy as np

from collections import defaultdict

hicFile = sys.argv[1]

refFile = sys.argv[2]
binSize = int(sys.argv[3])

print("#####",file=sys.stderr)
print('using bin size of: {}.'.format(sys.argv[3]),file=sys.stderr)
print('binning .hic file: {}'.format(sys.argv[1]),file=sys.stderr)
print('using reference digest file: {}'.format(sys.argv[2]),file=sys.stderr)
          
###first, load the reference
chrList = []

# Frags contains the fragment start positions, while fragPoss contains the fragment number for that fragment
frags = {}
fragToPoss = {}
convertFragToPos = {}

curPoss = 0
for entry in open(refFile):
    l = entry.strip().split()
    L = len(l[1:])
    chrList.append(l[0])
    frags[l[0]] = list(map(int,l[1:]))
    fragToPoss[l[0]] = list(range(curPoss, curPoss + L ))
    for pos in range(L):
        convertFragToPos[fragToPoss[l[0]][pos]] = (l[0], l[1+pos])
    curPoss += L

print("original fragment position summary")
for c, data in frags.items():
    print('{}: {} ... {}'.format(c, ' '.join(list(map(str,data[:3]))), ' '.join(list(map(str,data[-3:])))))

print("original fragment position index summary") #this table corresponds to the values in .hic entries that need to be mapped to a bin   
for c, data in fragToPoss.items():
    print('{}: {} ... {}'.format(c, ' '.join(list(map(str,data[:3]))), ' '.join(list(map(str,data[-3:])))))

#key on chr, values are start stop positions of bins
bins = defaultdict(list)
curBinNum = 0
#bin data must be constructed based on fragment positions
for ch, data in frags.items():
    curBin = binSize
    for pos in data:
        if pos > curBin:
            curBin += binSize
            curBinNum += 1
        bins[ch].append(curBinNum)
    curBinNum += 1

chBoundaries = []
curBin = 0
for c,data in bins.items():
    curBin = data[-1] + 1
    #the -0.5 pushes the dashes to between the bins as they should be
    chBoundaries.append(curBin-0.5)
    

print(chBoundaries)
###compute the mapping of fragment to bin into a dict

print("Bin level position summary")        
for c, data in bins.items():
    print('{}: {} ... {}'.format(c, ' '.join(list(map(str,data[:3]))), ' '.join(list(map(str,data[-3:])))))

binToPosition = {}
for c, data in bins.items():
    L = len(data)
    start = 0
    pos = 0
    while pos < L-1:
        while pos < L -1 and data[pos] == data[pos+1]:
            pos += 1
        if pos < L -1:
            end = frags[c][pos]
            #print( data[pos], (c,start,end))
            binToPosition[data[pos]] = (c,start,end)
            start = end
            pos += 1
        else:
            binToPosition[data[pos]] = (c,start,frags[c][-1])

#this may need to be changed so that the bins are exact multiples of the binSize,
#   which is dodgy when the fragments are not multiplies of the bin sizes. A
#   rule must be asserted to determine how that discontinuity should be reconciled.

###calculate bin number -> position and save to "fragmap" file
refBinsOut = open(sys.argv[1].replace('.hic.txt','.bed'),'w')
for binNum, span in binToPosition.items():
    l1,l2,l3 = span
    refBinsOut.write('{}\t{}\t{}\t{}'.format(l1,l2,l3,binNum))
    
BinToPos = {}
fragToBin = {}

for ch in bins.keys():
    for pos in range(len(fragToPoss[ch])):
        fragToBin[fragToPoss[ch][pos]] = bins[ch][pos]
        BinToPos[bins[ch][pos]] = (ch, fragToPoss[ch][pos])
        
totalBinCount = curBinNum
        
#for fragID, bin_ in fragToBin.items():
#    print( '{}:{}'.format(fragID,bin_))
###load the hic file, storing the bin number for each fragment pair (converted using results from above)
hicData = []
for entry in open(hicFile):
    l = entry.strip().split()
    #convert all but the readID into int type
    for pos in [1,3,4,5,7,8,9,10]:
        l[pos] = int(l[pos])
    hicData.append(l)



contactMatrix = np.zeros((curBinNum, curBinNum), dtype='u2')
#going below about 20kb resolution will require a sparse matrix solution

#convert hic fragments into bins and populate the results into a matrix
for entry in hicData:
    x = fragToBin[entry[4]] 
    y = fragToBin[entry[8]]
    contactMatrix[x,y] += 1
    contactMatrix[y,x] += 1

matrixOut = open(sys.argv[1].replace('.hic.txt','.matrix'),'w')

for x in range(curBinNum-1):
    for y in range(x, curBinNum):
        if contactMatrix[y,x] > 0:
            matrixOut.write('{}\t{}\t{}'.format(y,x,contactMatrix[y,x]))
