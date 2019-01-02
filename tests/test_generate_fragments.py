import pytest
import re
from pore_c.tools.generate_fragments import create_regex, revcomp, find_site_positions


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

