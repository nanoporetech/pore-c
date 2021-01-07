import re
from pathlib import Path
from typing import List, Pattern, Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from pore_c.datasources import IndexedFasta
from pore_c.utils import kmg_bases_to_int

from ..config import PQ_ENGINE, PQ_VERSION
from ..model import FragmentRecord, FragmentRecordDf


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
    reference_fasta: Path, digest_type: str, digest_param: str, fragments: Path = None, digest_stats: Path = None
) -> Tuple[FragmentRecordDf, pd.DataFrame]:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""

    # convert the sequences to a dask bag
    ff = IndexedFasta(reference_fasta)
    seq_bag = ff.to_dask()
    dtype = FragmentRecord.pandas_dtype(overrides={"chrom": pd.CategoricalDtype(ff._chroms, ordered=True)})

    frag_df = (
        pd.concat(
            seq_bag.map(lambda x: (x["seqid"], x["seq"], digest_type, digest_param))
            .starmap(create_fragment_dataframe)
            .compute()
        )
        .astype(dict(chrom=dtype["chrom"]))
        .sort_values(["chrom", "start"])
        .assign(fragment_id=lambda x: np.arange(len(x), dtype=int) + 1)
        .astype(dtype)
    )

    frag_df.to_parquet(str(fragments), engine=PQ_ENGINE, index=False, version=PQ_VERSION)

    summary_stats_df = (
        frag_df.groupby("chrom")["fragment_length"]
        .agg(["size", "mean", "median", "min", "max"])
        .fillna(-1)
        .astype({"size": int, "min": int, "max": int})
        .rename(columns={"size": "num_fragments"})
    )
    summary_stats_df.to_csv(digest_stats)

    return frag_df, summary_stats_df


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
    except Exception as exc:
        raise ValueError(
            "Error compiling regex for pattern {}, redundance form: {}\n{}".format(pattern, fwd_rev_pattern, exc)
        )
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
        positions = list(range(0, len(seq), bin_width))
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
