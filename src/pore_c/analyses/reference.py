import re
from pathlib import Path
from time import sleep
from typing import Iterator, List, NamedTuple, Pattern

import dask.dataframe as dd
import numpy as np
import pandas as pd
import yaml
from intake import open_catalog
from intake.catalog.local import YAMLFileCatalog
from pandas import DataFrame
from pysam import FastaFile

from pore_c.datasources import IndexedFasta
from pore_c.model import FragmentDf
from pore_c.utils import kmg_bases_to_int

# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans("ACGTWSMKRYBDHVNacgtwsmkrybdhvn-)(][", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn-()[]")

# translate degenerate bases to regex set
DEGENERATE_TRANS = str.maketrans(
    dict(
        [
            ("N", "[ACGT]"),
            ("V", "[ACG]"),
            ("H", "[ACT]"),
            ("D", "[AGT]"),
            ("B", "[CGT]"),
            ("W", "[AT]"),
            ("S", "[CG]"),
            ("M", "[AC]"),
            ("K", "[GT]"),
            ("R", "[AG]"),
            ("Y", "[CT]"),
        ]
    )
)


def create_virtual_digest(
    reference_fasta: IndexedFasta,
    digest_type: str,
    digest_param: str,
    fragment_df_path: Path,
    summary_stats_path: Path,
    n_workers: int = 1,
) -> FragmentDf:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""

    parallel = n_workers > 1
    if parallel:
        from dask.distributed import Client, LocalCluster

        cluster = LocalCluster(processes=True, n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    # convert the sequences to a dask bag
    seq_bag = reference_fasta.to_dask()
    chrom_dtype = pd.CategoricalDtype(reference_fasta._chroms, ordered=True)
    FragmentDf.set_dtype("chrom", chrom_dtype)

    frag_df = (
        pd.concat(
            seq_bag.map(lambda x: (x["seqid"], x["seq"], digest_type, digest_param))
            .starmap(create_fragment_dataframe)
            .compute()
        )
        .astype({"chrom": chrom_dtype})
        .sort_values(["chrom", "start"])
        .assign(fragment_id=lambda x: np.arange(len(x), dtype=int) + 1)
        .fragmentdf.cast(subset=True)
    )

    if parallel:
        while True:
            processing = client.processing()
            still_running = [len(v) > 0 for k, v in processing.items()]
            if any(still_running):
                sleep(10)
            else:
                break
        client.close()
        cluster.close()

    # use pandas accessor extension
    frag_df.fragmentdf.assert_valid()

    frag_df.to_parquet(str(fragment_df_path), index=False)

    summary_stats = (
        frag_df.groupby("chrom")["fragment_length"]
        .agg(["size", "mean", "median", "min", "max"])
        .fillna(-1)
        .astype({"size": int, "min": int, "max": int})
        .rename(columns={"size": "num_fragments"})
    )
    summary_stats.to_csv(summary_stats_path)

    return frag_df


def revcomp(seq: str) -> str:
    """Return the reverse complement of a string:
    """
    return seq[::-1].translate(COMPLEMENT_TRANS)


def replace_degenerate(pattern: str) -> str:
    """Replace degenerate bases with regex set"""
    return pattern.translate(DEGENERATE_TRANS)


def create_regex(pattern: str) -> Pattern:
    """Takes a raw restriction digest site in the form of a regular expression
    string and returns a regular expression object consisting of both forward
    and reverse complement versions of the pattern"""

    site_raw = pattern.replace("(", " ").replace(")", " ").replace(" ", "").replace("|", " ").split()
    sites_raw = []
    for entry in sites_raw:
        sites_raw.append(revcomp(entry))

    sites_raw += site_raw
    if len(sites_raw) > 1:
        fwd_rev_pattern = "(" + "|".join(sorted(list(set(sites_raw)))) + ")"
    else:
        fwd_rev_pattern = "(" + sites_raw[0] + ")"

    ###
    fwd_rev_pattern = replace_degenerate(fwd_rev_pattern)

    try:
        regex = re.compile(fwd_rev_pattern)
    except:
        raise ValueError("Error compiling regex for pattern {}, redundance form: {}".format(pattern, fwd_rev_pattern))
    return regex


def find_fragment_intervals(digest_type: str, digest_param: str, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    if digest_type == "regex":
        regex = create_regex(digest_param)
        positions = find_site_positions_regex(regex, seq)
    elif digest_type == "bin":
        bin_width = kmg_bases_to_int(digest_param)
        positions = find_site_positions_bins(bin_width, seq)
    elif digest_type == "enzyme":
        positions = find_site_positions_biopython(digest_param, seq)
    intervals = to_intervals(positions, len(seq))
    return intervals


def to_intervals(positions: List[int], chrom_length: int):
    prefix, suffix = [], []
    if (len(positions) == 0) or positions[0] != 0:
        prefix = [0]
    if (len(positions) == 0) or positions[-1] != chrom_length:
        suffix = [chrom_length]
    endpoints = np.array(prefix + positions + suffix)
    return {"start": endpoints[:-1], "end": endpoints[1:]}


def find_site_positions_bins(bin_width, seq: str) -> List[int]:
    """Mimic a fixed-width sequence digest by returning the positions of fixed-width bin boundaries"""
    if len(seq) < bin_width:
        return []
    else:
        positions = list(range(bin_width, len(seq), bin_width))
        return positions


def find_site_positions_regex(regex: Pattern, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    positions = [m.start() for m in regex.finditer(seq.upper())]
    return positions


def find_site_positions_biopython(enzyme: str, seq: str) -> List[int]:
    from Bio import Restriction
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        raise ValueError("Enzyme not found: {}".format(enzyme))
    s = Seq(seq, IUPACAmbiguousDNA())
    positions = [_ - 1 for _ in enz.search(s)]
    return positions


def create_fragment_dataframe(seqid: str, seq: str, digest_type: str, digest_param: str) -> DataFrame:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""
    intervals = (
        DataFrame(find_fragment_intervals(digest_type, digest_param, seq))
        .assign(chrom=seqid)
        .eval("fragment_length = end - start")
    )
    return intervals


def create_virtual_digest_dataframe(reference_fasta: str, digest_type: str, digest_param) -> pd.DataFrame:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""

    return fragment_df
