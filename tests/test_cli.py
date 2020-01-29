import os
from pathlib import Path

from click.testing import CliRunner

from pore_c.cli import cli


# def test_logging():
#    runner = CliRunner(mix_stderr=False)
#    result = runner.invoke(cli, ['-vvv', 'refgenome', 'catalog'])
#    assert 'pore_c - DEBUG - Logger set up' in result.stdout


def test_prepare_refgenome(raw_refgenome_file):
    runner = CliRunner()

    with runner.isolated_filesystem():
        cwd = Path(os.getcwd())
        result = runner.invoke(cli, ["refgenome", "prepare", str(raw_refgenome_file), "refgenome"])
        assert result.exit_code == 0
        new_files = set([f.name for f in cwd.glob("*.*")])
        expected_files = set(
            [
                "refgenome.chromsizes",
                "refgenome.metadata.csv",
                "refgenome.fa.fai",
                "refgenome.catalog.yaml",
                "refgenome.fa",
            ]
        )
        assert new_files == expected_files


def test_virtual_digest(raw_refgenome_file):
    runner = CliRunner()

    with runner.isolated_filesystem():

        result = runner.invoke(cli, ["refgenome", "prepare", str(raw_refgenome_file), "refgenome"])
        cwd = Path(os.getcwd())
        result = runner.invoke(cli, ["refgenome", "virtual-digest", "refgenome.fa", "NlaIII", "digest"])
        assert result.exit_code == 0
        new_files = set([f.name for f in cwd.glob("*.*")])
        expected_files = set(
            [
                "digest.catalog.yaml",
                "digest.digest_stats.csv",
                "digest.fragments.parquet",
                "refgenome.chromsizes",
                "refgenome.metadata.csv",
                "refgenome.fa.fai",
                "refgenome.catalog.yaml",
                "refgenome.fa",
            ]
        )
        assert new_files == expected_files
