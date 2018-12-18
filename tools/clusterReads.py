#!/usr/bin/env python

# a cheap and easy way to cluster mapping hits.
# assumes discrete mapping of monomers within X bases

import sys
import pysam as ps
import networkx as nx

bamIn = ps.AlignmentFile(sys.argv[1])
bamOut = ps.AlignmentFile(sys.argv[1].replace(".bam",".split.{}.bam".format(sys.argv[2])),'wb',template=bamIn)
removedOut = ps.AlignmentFile(sys.argv[1].replace(".bam",".split.{}.removed.bam".format(sys.argv[2])),'wb',template=bamIn)

trimSize = int(sys.argv[2])

__doc__ = """
Takes all reads from a namesorted bam file. trims back each read by trimSize 
number of bases, then clusters the reads by overlap. If reads overlap with 
each other they get placed in a cluster together. From each cluster, the highest 
scoring alignment is included in the filtered data set."""
mappings = []
tot_count = 0
inc_count = 0
exc_count = 0
lastRead = False

for entry in bamIn:
    tot_count += 1
    if not lastRead:
        lastRead = entry.query_name
        mappings = [entry]
    elif lastRead == entry.query_name:
        mappings.append(entry)
    elif lastRead != entry.query_name:
        if mappings[0].is_unmapped:
            removedOut.write(mappings[0])
            exc_count += 1
            lastRead = entry.query_name
            mappings = [entry]
            continue
        L = len(mappings)
        #reads with only one mapping are not helpful and are here culled
        if L > 1:
            G = nx.Graph()
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
                bestClustScore = 0
                bestRead = -1
                for node in clust:
                    score = dict(mappings[node].tags)["AS"]
                    if score > bestClustScore:
                        bestRead = node
                        bestClustScore = score
                bamOut.write(mappings[bestRead])
                inc_count += 1
                for node in clust:
                    if node != bestRead:
                        removedOut.write(mappings[node])
                        exc_count += 1
        else:#if single mapping
            removedOut.write(mappings[0])
            exc_count += 1
        mappings = [entry]
        lastRead = entry.query_name
### and one last graph partitioning for the last read cluster, which has no read after it
    
if not mappings[0].is_unmapped:
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
        if len(Z) > 1:
            for clust in Z:
                readsOut = 0
                bestClustScore = 0
                bestRead = -1
                for node in clust:
                    score = dict(mappings[node].tags)["AS"]
                    if score > bestClustScore:
                        bestRead = node
                        bestClustScore = score
                bamOut.write(mappings[bestRead])
                inc_count += 1
                for node in clust:
                    if node != bestRead:
                        removedOut.write(mappings[node])
                        exc_count += 1
        else:
            for node in Z[0]:
                removedOut.write(mappings[node])
    else:
        removedOut.write(mappings[0])
else:
    removedOut.write(mappings[0])

print("{} total. {} included. {} excluded.".format(tot_count,inc_count,exc_count),file = sys.stderr)
