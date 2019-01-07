import pytest
from pore_c.tools.map_to_frags import load_reference, assign_to_fragments, porec_iterator, map_to_fragments


def test_load_reference(hicREF_file):
    ref_sizes, ref_IDs = load_reference(hicREF_file)

    #the sizes dictionary and the IDs dictionary have the same number of keys and values
    assert ref_sizes.keys() == ref_IDs.keys()
    assert len(ref_sizes.values()) == len(ref_IDs.values())
    IDs_collected = []

    #all fragment IDs are unique
    for ch, IDs in ref_IDs.items():
        IDs.collected.extend(IDs)
    assert len(set(IDs_collected)) == len(IDs_collected)


    for positions in ref_sizes.values():
        assert sorted(positions) == positions
    assert len(set(IDs_collected)) == len(IDs_collected)


    #that the 

@pytest.mark.parametrize(" ", { })
def test_assign_to_fragments(hicREF_file):
    pass

def test_porec_iterator():
    pass

def test_map_to_fragments():
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
