import pytest
from pysam import AlignmentFile
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="function")
def namesorted_align_file():
    af = AlignmentFile(DATA_DIR / "test_ns.sam")
    return af

