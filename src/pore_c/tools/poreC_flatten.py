import bisect
import gzip
import re
import sys
import uuid
from collections import defaultdict
from itertools import combinations, combinations_with_replacement

import numpy as np
from pysam import AlignmentFile


class Contact:
    def __init__(self, ch, frag, strand, poss, mapq):
        self.ch = ch
        self.fragID = int(frag)
        if type(strand) == bool:
            self.strand = 16 if strand else 0
        else:
            self.strand = strand
        self.poss = poss
        self.mapq = mapq

    def __str__(self):
        return " ".join(
            list(map(str, (self.ch, self.fragID, self.strand, self.poss, self.mapq)))
        )


class Cwalk:
    def __init__(self, readname=None):
        self.name = readname
        #        self.fragIDs = set()
        self.contacts = []

    def from_entry(self, entry):
        l = entry.strip().split()
        L = int((len(l[1:])) / 5)
        self.name = l[0]
        mapqs = l[-L:]
        entries = []

        for x in range(1, 4 * L, 4):
            entries.append(l[x : x + 4])

        for x in range(L):
            strand, ch, poss, frag = entries[x]
            if strand == "16":
                strand = True
            else:
                strand = False
            self.add(Contact(ch, frag, strand, poss, mapqs[x]))

    def add(self, monomer):
        self.contacts.append(monomer)

    def sort(self):
        self.contacts.sort(key=lambda c: c.fragID)

    def length(self):
        return len(self.contacts)

    def get_chrs(self):
        return tuple([x.ch for x in self.contacts])

    # returns a dictionary of lists of cwalk objects keyed on sorted chromosome sets of size=size
    def flatten(self, size, direct=False):
        outputs = defaultdict(list)
        chrList = self.get_chrs()
        L = len(self.contacts)

        if len(self.contacts) >= size:
            if direct:
                for pos in range(len(self.contacts) - size):
                    contact_group = self.contacts[pos : pos + size]
                    new_walk = Cwalk("{}_{}".format(self.name, pos))
                    chrs = tuple(sorted([x.ch for x in contact_group]))
                    for contact in contact_group:
                        new_walk.add(contact)
                    outputs[chrs].append(new_walk)
                return outputs

            else:
                for c, contact_group in enumerate(combinations(self.contacts, r=size)):
                    new_walk = Cwalk("{}_{}".format(self.name, c))
                    chrs = tuple(sorted([x.ch for x in contact_group]))
                    for contact in contact_group:
                        new_walk.add(contact)
                    outputs[chrs].append(new_walk)
            return outputs
        else:
            return None

    #        c = 0
    #        if len(self.contacts) >= size:
    #            for contact_group in combinations(range(L),r = size):
    #                new_walk = Cwalk("{}_{}".format(self.name,c))
    #                chrs = tuple(sorted([chrList[x] for x in contact_group]))
    #                for pos in contact_group:
    #                    new_walk.add(self.contacts[pos])
    #                outputs[chrs].append(new_walk)
    #                c += 1
    #            return outputs
    #        else:
    #            return None
    ###
    def __str__(self):
        mapString = []
        mapqs = []
        for x in self.contacts:
            tempString = "{strand} {ch} {pos} {frag}".format(
                strand=x.strand, ch=x.ch, pos=x.poss, frag=x.fragID
            )
            mapString.append(tempString)
            mapqs.append(x.mapq)

        mapString = " ".join(mapString)
        quals = " ".join(list(map(str, mapqs)))
        return "{name} {mappings} {quals}".format(
            name=self.name, mappings=mapString, quals=quals
        )

#file    format  type    num_seqs        sum_len min_len avg_len max_len Q1      Q2      Q3      sum_gap N50
#calculates contact_N50, mean, median,  min, max, Q1, Q2, Q3, total_monomer_count, total_entries for a .poreC file
def stats(poreC_file_in: str) -> str:
    entries = []
    contacts = []
    for entry in open(poreC_file_in):
        c = Cwalk("None")
        c.from_entry(entry)
        z = c.length()
        entries.append(z)
        contacts.append((z**2 - z) / 2)

    #reverse sort for N50 calculation using a[::-1].sort()
    entries = np.array(entries,dtype = float)
    entries[::-1].sort()

    contacts = np.array(contacts,dtype = float)
    contacts[::-1].sort()

    template = "{fn},{num_seqs},{total_monomer_count},{min_len},{mean_len},{median_len},{mode_len},{max_len},{Q1},{Q2},{Q3},{monomer_n50},{total_contact_count},{c_min_len},{c_mean_len},{c_median_len},{c_mode_len},{c_max_len},{c_Q1},{c_Q2},{c_Q3},{contact_n50}"
    print(template.replace("{","").replace("}","")) #print a header

    num_seqs = len(entries)
    if num_seqs == 0:
        print( template.format(fn = poreC_file_in,
                               num_seqs = num_seqs,
                               total_monomer_count = 0,
                               min_len = 0,
                               mean_len = 0,
                               mode_len = 0,
                               median_len = 0,
                               max_len = 0,
                               Q1 = 0,
                               Q2 =0,
                               Q3 =0,
                               poreC_n50 = 0)
               )
        exit()

    total_monomer_count = entries.sum()
    max_len = entries.max()
    min_len = entries.min()
    mean_len = entries.mean()
    median_len = np.median(entries)
    half_sum = total_monomer_count / 2.0
    temp_sum = 0
    monomer_n50 = 0
    while temp_sum < half_sum:
        temp_sum += entries[monomer_n50]
        monomer_n50 += 1

    (_, idx, counts) = np.unique(entries, return_index=True, return_counts=True)
    index = idx[np.argmax(counts)]
    mode_len = entries[index]


    total_contact_count = contacts.sum()
    c_max_len = contacts.max()
    c_min_len = contacts.min()
    c_mean_len = contacts.mean()
    c_median_len = np.median(contacts)
    half_sum = total_contact_count / 2.0
    temp_sum = 0

    temp_sum = 0
    contact_n50 = 0
    while temp_sum < half_sum:
        temp_sum += contacts[contact_n50]
        contact_n50 += 1


    (_, idx, counts) = np.unique(contacts, return_index=True, return_counts=True)
    index = idx[np.argmax(counts)]
    c_mode_len = contacts[index]

    print( template.format(fn = poreC_file_in,
                           num_seqs = num_seqs,
                           total_monomer_count = total_monomer_count ,
                           min_len = min_len,
                           mean_len = mean_len,
                           mode_len = mode_len,
                           median_len = median_len,
                           max_len = max_len,
                           Q1 = np.percentile(entries,25),
                           Q2 = np.percentile(entries,50),
                           Q3 = np.percentile(entries,75),
                           monomer_n50 = entries[monomer_n50],
                           total_contact_count = total_contact_count ,
                           c_min_len = c_min_len,
                           c_mean_len = c_mean_len,
                           c_mode_len = c_mode_len,
                           c_median_len = c_median_len,
                           c_max_len = c_max_len,
                           c_Q1 = np.percentile(contacts,25),
                           c_Q2 = np.percentile(contacts,50),
                           c_Q3 = np.percentile(contacts,75),
                           contact_n50 = contacts[contact_n50]))
                                    



# creates output filehandle
# reads input file
# passes off operations to the functions above as necessary
def flatten_multiway(file_in, file_out, size, sort=True, direct=False):
    # gather data by chromosome combination,
    #   so that it can be printed out in sorted order
    #   in order to play nicely with juicer

    # this value stores a list of chromosomes visited
    #   so that they can be written out in an order
    #   that is most compatible with other tools in the
    #   community.
    chr_used = set()
    contact_bundle = defaultdict(list)
    f_out = open(file_out, "w")

    with open(file_in) as f:
        for entry in f:
            walk = Cwalk()
            walk.from_entry(entry)
            if walk.length() < size:
                continue
            if sort:
                walk.sort()
            for contact in walk.contacts:
                chr_used.add(contact.ch)
            flattened_contacts = walk.flatten(size, direct=direct)
            # this is better than a = {**a, **b} because |entries| <<< |chrs|
            for c, entries in flattened_contacts.items():
                # contacts must be sorted into blocks (e.g. all chr1 intra-chr, then chr1-chr2, chr1-chr3,etc.) which is different from a simple sort.
                for contact in entries:
                    f_out.write(str(contact) + "\n")

    f_out.close()


def make_salsa_bedfile(
    hictxt_file_in: str, bedfile_out: str, frag_bed_ref: str
) -> None:
    fragments = {}
    for entry in gzip.open(frag_bed_ref, "rt"):
        l = entry.strip().split()
        fragments[l[3]] = l

    fOut = open(bedfile_out, "w")

    entry_template = "{ch}\t{start}\t{end}\t{read_id}/{pairID}\t{mapq}\t{strand}\n"

    for entry in open(hictxt_file_in):
        l = entry.strip().split()
        frag = fragments[l[4]]
        fOut.write(
            entry_template.format(
                ch=frag[0],
                start=frag[1],
                end=frag[2],
                read_id=l[0],
                pairID=1,
                mapq=l[9],
                strand="-" if l[1] else "+",
            )
        )
        frag = fragments[l[8]]
        fOut.write(
            entry_template.format(
                ch=frag[0],
                start=frag[1],
                end=frag[2],
                read_id=l[0],
                pairID=2,
                mapq=l[10],
                strand="-" if l[5] else "+",
            )
        )

    fOut.close()


#this code was moved to the stats method of the map-to-frags tool.
#There may be some use for a stand-alone tool at some point?
deprecated = """
def per_read_intra_inter(porec_file_in: str, csv_out: str, frag_bed_ref: str) -> None:
    entry_template = "{read_id},{intra},{inter},{ratio},{pct_intra}\n"
    header = "read_id,intra,inter,ratio,pct_intra\n"

    f_out = open(csv_out, "w")

    f_out.write(header)

    for entry in map(Cwalk.from_entry, open(porec_file_in)):
        intra = 0
        inter = 0
        for c1 in range(len(entry) - 1):
            for c2 in range(c1 + 1, len(entry)):
                if entry.contacts[c1].ch == entry.contacts[c2].ch:
                    intra += 1
                else:
                    inter += 1
        f_out.write(
            entry_template.format(
                read_id=entry.name,
                intra=intra,
                inter=inter,
                ratio=intra / float(inter),
                pct_intra=(intra / float(intra + inter)),
            )
        )
"""


def fragment_end_metrics(bam_file_in: str, csv_out: str, hicref: str):
    entry_template = (
        "{read_id},{frag_num},{startpoint_d_to_motif},{endpoint_d_to_motif}\n"
    )
    header = "read_id,frag_num,startpoint_d_to_motif,endpoint_d_to_motif\n"

    f_out = open(csv_out, "w")

    f_out.write(header)

    cutmap = {}
    for entry in open(hicref):
        l = entry.strip().split()
        cutmap[l[0]] = np.array(list(map(int, [0] + l[1:])), dtype=int)

    last_readname = False
    for entry in AlignmentFile(bam_file_in):

        start = bisect.bisect_left(cutmap[entry.reference_name], entry.reference_start)
        l_start = max(0, start - 1)
        r_start = min(len(cutmap[entry.reference_name]) - 1, start + 1)
        if len(cutmap[entry.reference_name]) > 3:
            d_start = min(
                map(
                    lambda x: abs(
                        cutmap[entry.reference_name][x] - entry.reference_start
                    ),
                    [l_start, start, r_start],
                )
            )
        else:
            d_start = min(
                map(
                    lambda x: abs(
                        cutmap[entry.reference_name][x] - entry.reference_start
                    ),
                    [0, 1],
                )
            )
        end = bisect.bisect_left(cutmap[entry.reference_name], entry.reference_end)
        l_end = max(0, end - 1)
        end = min(end, len(cutmap[entry.reference_name]) - 1)
        r_end = min(len(cutmap[entry.reference_name]) - 1, end + 1)
        if len(cutmap[entry.reference_name]) > 2:
            try:
                d_end = min(
                    map(
                        lambda x: abs(
                            cutmap[entry.reference_name][x] - entry.reference_end
                        ),
                        [l_end, end, r_end],
                    )
                )
            except:
                print("oops:", [l_end, end, r_end], cutmap[entry.reference_name][x])
        else:
            d_start = min(
                map(
                    lambda x: abs(
                        cutmap[entry.reference_name][x] - entry.reference_end
                    ),
                    [0, 1],
                )
            )

        if not last_readname or last_readname != entry.query_name:
            idx = 0
        else:
            idx += 1

        f_out.write(
            entry_template.format(
                read_id=entry.query_name,
                frag_num=idx,
                startpoint_d_to_motif=d_start,
                endpoint_d_to_motif=d_end,
            )
        )

        last_readname = entry.query_name


