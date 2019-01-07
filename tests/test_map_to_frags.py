import pytest
from pore_c.tools.map_to_frags import load_reference, assign_to_fragments, porec_iterator, map_to_fragments
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
def test_assign_to_fragments(hicREF_file, namesorted_align_file):
    ref_sizes, ref_IDs = load_reference(hicREF_file)

    results = []
    for readID, aligns in porec_iterator(namesorted_align_file):
        walk = Cwalk(name)
        for monomer in aligns:
            loc = (monomer.reference_name, monomer.reference_start,monomer.reference_end)
            frag,pos = assign_to_fragment(ref_frags, ref_IDs, loc, method = method)
            new_contact = Contact(monomer.reference_name, frag, monomer.is_reverse, pos, monomer.mapping_quality)
            walk.add(new_contact)
        results.append(len(str(walk).split()))

    assert results == [11, 26, 21, 21, 11, 6]

def test_porec_iterator(namesorted_align_file):
    res = dict([aligns for aligns in porec_iterator(namesorted_align_file)])
    assert res.keys() == ["0a4cf95d-1410-42b1-867c-68e987123675","0a6b6929-f064-42d4-afd4-24b38a966d03","0a0030d2-60a6-4498-8ea0-329d9dcac82b","0aa14eb3-bdd8-40eb-97c3-d37c6eb99ef6","0af7db4f-2804-4dc8-aae9-fc4611008d4a","0b72d767-a437-43e2-a9ea-717d99864688"]
    assert([len(align) for align in res.values()] = [2, 5, 6, 5, 2, 1])

    #check that all the alignments in each read bundle have the same query_name
    # and that the order in query is preserved (i.e. that they're sorted by query_start)
    for readID, aligns in res.items():
        q_name = aligns[0].query_name
        #bundle names are all identical
        assert( all([x.query_name == q_name for x in aligns]))

    for aligns in res.values():
        #sorted order
        assert( sorted(aligns,key=lambda x: x.query_alignment_start) == aligns)

#signature: def map_to_fragments(input_bam: str, ref_file: str, output_file: str, method: str) -> None:
def test_map_to_fragments(namesorted_align_file, hicREF_file, 'start'):
    pass


#def test_read_mappings_iter(namesorted_align_file):
#    res = [len(aligns) for aligns in read_mappings_iter(namesorted_align_file)]
#    assert(res == [2, 5, 6, 5, 2, 1])


#def test_cluster_aligned_segments(namesorted_align_file):

#    aligns = list(namesorted_align_file)[2:7]
#    keep = cluster_aligned_segments(aligns, 20)
#    assert(keep == [0, 1, 3])


@pytest.mark.parametrize("pattern,rev_pattern",[
    ('ACGT', 'ACGT'),
    ('ACAGTA', 'TACTGT'),
    ('(ACGT)', '(ACGT)'),
    ('(A-CGT)', '(ACG-T)'),
    ('AAGCTT', 'AAGCTT'),
])
def test_revcomp(pattern, rev_pattern):
    assert(rev_pattern == revcomp(pattern))


@pytest.mark.parametrize("pattern,regex",[
    ('AAGCTT', '(AAGCTT|AAGCTT)'),
    ('AAGNCTT', '(AAG[ACGT]CTT|AAG[ACGT]CTT)'),
    ('(GAATTC|GCGGCCGC)', '((GAATTC|GCGGCCGC)|(GCGGCCGC|GAATTC))')
])
def test_create_regex(pattern, regex):
    assert(regex == create_regex(pattern).pattern)


@pytest.mark.parametrize("regex,seq,positions",[
    ('(ACAGTA|TACTGT)', "ACAGTATACTGT", [0, 6]),
    ('(ACAGTA|TACTGT)', "ACAGAATACAGT", []),
    ('(AGCT|AAGCTT)', 'AAGCTT', [0]),
    ('(AAGCTT|AAGCTT)', 'AAGCTT', [0]),
    ('(AAGCTT|AAGCTT)', 'NAAGCTT', [1]),
    ('(AAGCTT|AAGCTT)', 'AAGCTTAAGCTT', [0,6]),
])
def test_find_positions(regex, seq, positions):
    assert(positions == find_site_positions(re.compile(regex), seq))
