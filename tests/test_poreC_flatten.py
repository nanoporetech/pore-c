import pytest

from pore_c.tools.poreC_flatten import flatten_multiway as flatten_multiway
from pore_c.tools.map_to_frags import Contact, Cwalk

@pytest.mark.parametrize("poreC,size,sort,monomers",[ ("00aaefda-bc3e-4e34-9cab-3ef38f59c5eb 0 1 160205974 38790 16 1 182085591 45470 16 1 192474669 48773 16 1 192480456 48775 42 248 44 26",2,True,
['00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_0 0 1 160205974 38790 16 1 182085591 45470 42 248',
'00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_1 0 1 160205974 38790 16 1 192474669 48773 42 44',
'00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_2 0 1 160205974 38790 16 1 192480456 48775 42 26',
'00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_3 16 1 182085591 45470 16 1 192474669 48773 248 44',
'00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_4 16 1 182085591 45470 16 1 192480456 48775 248 26',
 '00aaefda-bc3e-4e34-9cab-3ef38f59c5eb_5 16 1 192474669 48773 16 1 192480456 48775 44 26'])
])
def test_Cwalk_flatten(poreC,size, sort,monomers):
    c = Cwalk()
    c.from_entry(poreC)
    if sort:
        c.sort()
    results = []
    pre_results = c.flatten(size)
    for entries in pre_results.values():
        results.extend([str(_).strip() for _ in entries])
    assert monomers == results

