import pytest
from pore_c.tools.cluster_reads import read_mappings_iter, cluster_aligned_segments


def test_read_mappings_iter(namesorted_align_file):
    res = [len(aligns) for aligns in read_mappings_iter(namesorted_align_file)]
    assert(res == [2, 5, 6, 5, 2, 1])


def test_cluster_aligned_segments(namesorted_align_file):

    aligns = list(namesorted_align_file)[2:7]
    keep = cluster_aligned_segments(aligns, 20)
    assert(keep == [0, 1, 3])


