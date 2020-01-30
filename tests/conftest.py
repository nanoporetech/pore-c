from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def raw_refgenome_file():
    return DATA_DIR / "GRCh38.fasta.gz"


@pytest.fixture(scope="session")
def read_fastq_file():
    return DATA_DIR / "GM12878_NlaIII_reads.fq.gz"


@pytest.fixture(scope="session")
def read_fastq_table():
    return DATA_DIR / "GM12878_NlaIII_reads.tsv"


@pytest.fixture(scope="session")
def read_sorted_bam():
    return DATA_DIR / "NlaIII_run01_GRCh38.read_sort.bam"
