import re
from pathlib import Path
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
from pore_c.model import FragmentMap, SeqDigest
from pore_c.utils import kmg_bases_to_int

# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn-)(][", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn-()[]"
)

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

    site_raw = (
        pattern.replace("(", " ").replace(")", " ").replace(" ", "").replace("|", " ").split()
    )
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
        raise ValueError(
            "Error compiling regex for pattern {}, redundance form: {}".format(
                pattern, fwd_rev_pattern
            )
        )
    return regex


def find_site_positions(regex_or_enzyme: str, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence"""
    if regex_or_enzyme.startswith("regex:"):
        regex = create_regex(regex_or_enzyme.split("regex:", 1)[1])
        positions = find_site_positions_regex(regex, seq)
    elif regex_or_enzyme.startswith("bin:"):
        bin_width = kmg_bases_to_int(regex_or_enzyme.split("bin:", 1)[1])
        positions = find_site_positions_bins(bin_width, seq)
    else:
        positions = find_site_positions_biopython(regex_or_enzyme, seq)
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


def create_fragment_dataframe(seqid: str, seq: str, restriction_pattern: str) -> DataFrame:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""
    intervals = (
        DataFrame(find_site_positions(restriction_pattern, seq))
        .assign(chrom=seqid)
        .eval("fragment_length = end - start")
    )
    return intervals


def create_virtual_digest_dataframe(
    reference_fasta: str, restriction_pattern: str = None
) -> pd.DataFrame:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""

    seq_bag = reference_fasta.to_dask()
    chrom_dtype = pd.CategoricalDtype(reference_fasta._chroms, ordered=True)

    fragment_df = (
        pd.concat(
            seq_bag.map(lambda x: (x["seqid"], x["seq"], restriction_pattern))
            .starmap(create_fragment_dataframe)
            .compute()
        )
        .astype({"chrom": chrom_dtype})
        .sort_values(["chrom", "start"])
        .assign(fragment_id=lambda x: np.arange(len(x), dtype=int))
        .loc[:, ["chrom", "start", "end", "fragment_id", "fragment_length"]]
    )
    return fragment_df


def create_virtual_digest(
    reference_fasta: IndexedFasta, restriction_pattern: str, output_prefix: Path
) -> pd.DataFrame:
    """Iterate over the sequences in a fasta file and find the match positions for the restriction fragment"""

    catalog_path = output_prefix.with_suffix(".virtual_digest.catalog.yaml")
    parquet_path = output_prefix.with_suffix(".virtual_digest.parquet")
    summary_path = output_prefix.with_suffix(".virtual_digest.summary.csv")

    for path in [catalog_path, parquet_path, summary_path]:
        if path.exists():
            raise IOError("Output path exists, please delete: {}".format(path))

    fragment_df = dd.from_pandas(
        create_virtual_digest_dataframe(reference_fasta, restriction_pattern), npartitions=1
    )
    fragment_df.to_parquet(str(parquet_path))

    summary_stats = (
        fragment_df.compute()
        .groupby("chrom")["fragment_length"]
        .agg(["size", "mean", "median", "min", "max"])
        .fillna(-1)
        .astype({"size": int, "min": int, "max": int})
        .rename(columns={"size": "num_fragments"})
    )
    summary_stats.to_csv(summary_path)

    catalog_data = {
        "name": "pore_c_virtual_digest",
        "description": "Output files of a pore-c tools virtual digest",
        "sources": {
            "fragment_df": {
                "driver": "parquet",
                "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(parquet_path.name)},
            },
            "summary_stats": {
                "driver": "csv",
                "args": {"urlpath": "{{ CATALOG_DIR }}/" + str(summary_path.name)},
            },
        },
    }
    with catalog_path.open("w") as fh:
        fh.write(yaml.dump(catalog_data))
    return YAMLFileCatalog(str(catalog_path))
