import pytest
import re
from pore_c.tools.map_to_frags import load_reference, assign_to_fragment, porec_iterator, map_to_fragments
from pore_c.tools.map_to_frags import Contact, Cwalk

def test_load_reference(hicREF_file):
    ref_sizes, ref_IDs = load_reference(hicREF_file)

    #the sizes dictionary and the IDs dictionary have the same number of keys and values
    assert ref_sizes.keys() == ref_IDs.keys()
    assert len(ref_sizes.values()) == len(ref_IDs.values())
    IDs_collected = []

    #all fragment IDs are unique
    for ch, IDs in ref_IDs.items():
        IDs_collected.extend(IDs)

    assert len(set(IDs_collected)) == len(IDs_collected)

    #that the position lists are loaded in sorted order
    for positions in ref_sizes.values():
        assert sorted(positions) == positions

#signature: assign_to_fragment(ref_frags: tuple, ref_IDs: dict, loc: tuple, method: str) -> NamedTuple:
#This is the function that actually does the assignment. map_to_fragments is a wrapper for it 
#  that manages bam reading and fileIO, whose internal code is here for management purposes

def test_assign_to_fragment(hicREF_file, namesorted_align_filename):
    ref_sizes, ref_IDs = load_reference(hicREF_file)

    results = []
    for readID, aligns in porec_iterator(namesorted_align_filename):
        walk = Cwalk(readID)
        for monomer in aligns:
            loc = (monomer.reference_name, monomer.reference_start,monomer.reference_end)
            frag,pos = assign_to_fragment(ref_sizes, ref_IDs, loc, 'start') 
            new_contact = Contact(monomer.reference_name, frag, monomer.is_reverse, pos, monomer.mapping_quality)
            walk.add(new_contact)
        results.append(len(str(walk).split()))

    assert results == [11, 26, 21, 21, 11, 6]


__doc__= """

(all positions indicated by brackets are inclusive)

                1          2         3
ref1: 123456789 0123456789 01234567890
      .........x..........x...........
   1:           [           ]           10,21       mapping overflows onto subsequent fragment
   2: [                  ]               1,19       mapping to two adjacent fragments

   4:                          [     ]  24,30       end of last fragment
   5:              [    ]               13,17       middle bit of a fragment
   6:           [        ]              10,19       perfect mapping 
                                 

              1         2     
ref2: 1234567 8901234567890 123 45
      .......x.............x...x..

   3: [     ]                            1, 6       beginning of first fragment,second chromosome


	

{"ref1" : [9,19,30], "ref2" : [7,20,23,25]}
{"ref1" : [0,1,2], "ref2" : [3,4,5,6]}
"""

#signature: assign_to_fragment(ref_frags: tuple, ref_IDs: dict, loc: tuple, method: str) -> NamedTuple:
@pytest.mark.parametrize("loc,frag_pos,frag_IDs,assign_method,fragID",[
(('ref1',10,21), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]}, {"ref1" : [0,1,2], "ref2" : [3,4,5,6]},"start",1),
(('ref1',1,19), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]},{"ref1" : [0,1,2], "ref2" : [3,4,5,6]},"start",0),
(('ref2',1,6), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]},{"ref1" : [0,1,2], "ref2" : [3,4,5,6]},"start",3),
(('ref1',24,30), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]},{"ref1" : [0,1,2], "ref2" : [3,4,5,6]},"start",2),
(('ref1',13,17), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]},{"ref1" : [0,1,2], "ref2" : [3,4,5,6]},"start",1),
(('ref1',10,19), {"ref1" : [9,19,30], "ref2" : [7,20,23,25]},{"ref1" : [0,1,2], "ref2" : [3,4,5,6]}, "start",1)]
)
def test_assign_to_fragment2(frag_pos,frag_IDs,loc,assign_method,fragID):
    ID,pos = assign_to_fragment(frag_pos,frag_IDs,loc,assign_method)
    assert ID == fragID



def test_porec_iterator(namesorted_align_filename):
    res = dict([aligns for aligns in porec_iterator(namesorted_align_filename)])
    assert sorted(res.keys()) == ['0a0030d2-60a6-4498-8ea0-329d9dcac82b', '0a4cf95d-1410-42b1-867c-68e987123675', '0a6b6929-f064-42d4-afd4-24b38a966d03', '0aa14eb3-bdd8-40eb-97c3-d37c6eb99ef6', '0af7db4f-2804-4dc8-aae9-fc4611008d4a', '0b72d767-a437-43e2-a9ea-717d99864688']
    assert([len(align) for align in res.values()] == [2, 5, 6, 5, 2, 1])

    #check that all the alignments in each read bundle have the same query_name
    # and that the order in query is preserved (i.e. that they're sorted by query_start)
    for readID, aligns in res.items():
        q_name = aligns[0].query_name
        #bundle names are all identical
        assert( all([x.query_name == q_name for x in aligns]))

    for aligns in res.values():
        #sorted order
        assert( sorted(aligns,key=lambda x: x.query_alignment_start) == aligns)

