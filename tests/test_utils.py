import pytest
from pore_c.utils import kmg_bases_to_int


@pytest.mark.parametrize("kmg_str,val",[
    ("5k", 5000),
    ("5m", 5000000),
    ("5g", 5000000000),
])
def test_kmgbases_to_int(kmg_str, val):
    assert(kmg_bases_to_int(kmg_str) == val)
