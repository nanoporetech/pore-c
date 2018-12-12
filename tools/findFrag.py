#!/usr/bin/env python


#takes in a location in a genome and a fragment file
# simply returns all fragments overlapping the location

import sys

frags = []

for entry in open(sys.argv[1]):
    l = entry.strip().split()
    frags.append((l[0],list(map(int,l[1:]))))


ch,pos = sys.argv[2].replace(',', '').split(':')
start,stop = pos.split('-')
start = int(start)
stop = int(stop)

pos = 0
ch_counter = 0
while ch != frags[ch_counter][0]:
    pos += len(frags[ch_counter][1])
    ch_counter += 1

ch_map = frags[ch_counter][1]
local_pos = 0
while ch_map[local_pos] < start:
    local_pos += 1

results = []

while ch_map[local_pos] < stop:
    results.append(pos + local_pos)
    local_pos += 1

for entry in results:
    print(entry)
