#!/usr/bin/env python

# a cheap and easy way to cluster mapping hits.
# assumes discrete mapping of monomers within X bases

import sys
import pysam as ps
import networkx as nx

bamIn = ps.AlignmentFile(sys.argv[1])
bamOut = ps.AlignmentFile(sys.argv[1].replace(".bam",".split.{}.bam".format(sys.argv[2])),'wb',template=bamIn)
removedOut = ps.AlignmentFile(sys.argv[1].replace(".bam","split.{}.removed.bam".format(sys.argv[2])),'wb',template=bamIn)

trimSize = int(sys.argv[2])

__doc__ = """
Takes all reads from a namesorted bam file. trims back each read by trimSize 
number of bases, then clusters the reads by overlap. If reads overlap with 
each other they get placed in a cluster together. From each cluster, the highest 
scoring alignment is taken into the Hi-C flattening step. This tool may or may
not be necessary based on the nature of the bwa alignments given the MC-4C 
parameter set. This needs to be confirmed before deprecating clusterReads.py
from the tool set."""

unmappedReads = 0
lastRead = False
for entry in bamIn:
    if entry.flag == 4:
        removedOut.write(entry)
        unmappedReads += 1
        continue
    if not lastRead:
        lastRead = entry.query_name
        G = nx.Graph()
        mappings = [entry]
    elif lastRead == entry.query_name:
        mappings.append(entry)
    else:
        print(entry.query_name)
        L = len(mappings)
        #reads with only one mapping are not helpful and are here culled
        if L > 1:
            for x in range(L - 1):
                G.add_node(x)
                for y in range(x+1, L):
                    #fudge the margins to prevent end to end clustering
                    x1 = min([mappings[x].query_alignment_start + trimSize, int((mappings[x].query_alignment_start +mappings[x].query_alignment_end) / 2 )])
                    x2 = max([mappings[x].query_alignment_end - trimSize,   int((mappings[x].query_alignment_start +mappings[x].query_alignment_end) / 2 )])
                    y1 = min([mappings[y].query_alignment_start + trimSize, int((mappings[y].query_alignment_start +mappings[y].query_alignment_end) / 2 )])
                    y2 = max([mappings[y].query_alignment_end - trimSize,   int((mappings[y].query_alignment_start +mappings[y].query_alignment_end) / 2 )])
                    if (x1 <= y1 and x2 <= y2 and x2 >= y1) or  (x1 >= y1 and x2 >= y2 and x2 <= y1):
                        #both dovetails
                        G.add_edge(x,y)
                    if (x1 <= y1 and x2 >= y2) or (x1 >= y1 and y2 >= x2):
                        #both cases of entirely contained
                        G.add_edge(x,y)
            G.add_node(L-1)
            Z = list(nx.connected_components(G))
            #a single cluster does not contain HiC data
            for clust in list(nx.connected_components(G)):
                readsOut = 0
                bestClustScore = 0
                bestRead = -1
                for node in clust:
                    if mappings[node].is_unmapped:
                        continue
                    score = dict(mappings[node].tags)["AS"]
                    if score > bestClustScore:
                        bestRead = node
                        bestClustScore = score
                if bestRead != -1:
                    #print(clust,bestRead)
                    bamOut.write(mappings[bestRead])
                    readsOut += 1
                    lastClust = bestRead
#                print ("Picking {} out of {} reads.".format(1, len(clust)))
            else:#if single mapping with no supplemental
                removedOut.write(mappings[0])
        mappings = [entry]
        lastRead = entry.query_name
        G = nx.Graph()

#...And one last graph partitioning for the last read cluster, 
#   which has no read after it to trigger evaluation while streaming through the file
        
L = len(mappings)
#reads with only one mapping are not helpful and are here culled
if L > 1:
    for x in range(L - 1):
        G.add_node(x)
        for y in range(x+1, L):
            #fudge the margins to prevent end to end clustering
            x1 = min([mappings[x].query_alignment_start + trimSize, int((mappings[x].query_alignment_start +mappings[x].query_alignment_end) / 2 )])
            x2 = max([mappings[x].query_alignment_end - trimSize,   int((mappings[x].query_alignment_start +mappings[x].query_alignment_end) / 2 )])
            y1 = min([mappings[y].query_alignment_start + trimSize, int((mappings[y].query_alignment_start +mappings[y].query_alignment_end) / 2 )])
            y2 = max([mappings[y].query_alignment_end - trimSize,   int((mappings[y].query_alignment_start +mappings[y].query_alignment_end) / 2 )])
            if (x1 <= y1 and x2 <= y2 and x2 >= y1) or  (x1 >= y1 and x2 >= y2 and x2 <= y1):
                #both dovetails
                G.add_edge(x,y)
            if (x1 <= y1 and x2 >= y2) or (x1 >= y1 and y2 >= x2):
                #both cases of entirely contained
                G.add_edge(x,y)
    G.add_node(L-1)
    Z = list(nx.connected_components(G))
    if len(Z) > 1:
        for clust in Z:
            readsOut = 0
            bestClustScore = 0
            bestRead = -1
            for node in clust:
                if mappings[node].is_unmapped:
                    removedOut.write(mappings[node])
                    continue
                score = dict(mappings[node].tags)["AS"]
                if score > bestClustScore:
                    bestRead = node
                    bestClustScore = score
            if bestRead != -1:
                bamOut.write(mappings[bestRead])
                readsOut += 1
                lastClust = bestRead
            rejects = set(clust).remove(bestRead)
            for rej in rejects:
                removedOut.write(mappings[rej])

    else:
        removedOut.write(mappings[Z[0]])
else:
    removedOut.write(mappings[0])
