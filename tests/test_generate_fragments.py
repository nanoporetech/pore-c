import pytest
from pore_c.tools.generate_fragments import create_regex, revcomp, find_site_positions, find_site_positions_biopython, find_site_positions_regex

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
    ('AAGCTT', '(AAGCTT)'),
    ('AAGNCTT', '(AAG[ACGT]CTT)'),
    ('(GAATTC|GCGGCCGC)', '(GAATTC|GCGGCCGC)')
])
def test_create_regex(pattern, regex):
    assert(regex == create_regex(pattern).pattern)


@pytest.mark.parametrize("regex,seq,positions",[
    ('(ACAGTA|TACTGT)', "ACAGTATACTGT", [0, 6]),
    ('(ACAGTA|TACTGT)', "ACAGAATACAGT", []),
    ('(ACAG[ACGT]A|TACTGT)', "ACAGAATACAGT", [0]),
    ('(ACAGTA|TAC[ACGT]GT)', "ACAGAATACAGT", [6]),
    ('(AGCT|AAGCTT)', 'AAGCTT', [0]),
    ('(AAGCTT)', 'AAGCTT', [0]),
    ('(AAGCTT)', 'NAAGCTT', [1]),
    ('(AAGCTT)', 'AAGCTTAAGCTT', [0,6]),
])
def test_find_positions_regex(regex, seq, positions):
    import re
    assert(positions == find_site_positions_regex(re.compile(regex), seq))


# HindIII cut positions (^ top strand, _ botton strand): A^AGCT_T
@pytest.mark.parametrize("enzyme,seq,positions",[
    ('HindIII', 'AAGCTT', [1]),
    ('HindIII', 'NAAGCTT', [2]),
    ('HindIII', 'AAGCTTAAGCTT', [1,7]),
])
def test_find_position_biopython(enzyme, seq, positions):
    assert(positions == find_site_positions_biopython(enzyme, seq))

# Note that these are the same site, but the regex gives the start of the pattern, not the cut site
@pytest.mark.parametrize("enzyme,seq,positions",[
    ('HindIII', 'AAGCTT', [1]),
    ('regex:AAGCTT', 'AAGCTT', [0]),
])
def test_find_position_facade(enzyme, seq, positions):
    assert(positions == find_site_positions(enzyme, seq))

