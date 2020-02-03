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


@pytest.fixture(scope="session")
def haplotagged_bam():
    return DATA_DIR / "NlaIII_run01_GRCh38.haplotagged.bam"


@pytest.fixture(scope="session")
def align_table_pq():
    return DATA_DIR / "NlaIII_run01_GRCh38.align_table.parquet"


@pytest.fixture(scope="session")
def fragment_table_pq():
    return DATA_DIR / "NlaIII_GRCh38.vd.fragments.parquet"
