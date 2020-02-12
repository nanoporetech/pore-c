from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def raw_refgenome_file():
    return DATA_DIR / "GRCh38.fasta.gz"
