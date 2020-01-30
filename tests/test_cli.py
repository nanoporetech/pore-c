import os
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from pore_c.cli import cli


def test_create_table(haplotagged_bam):
    runner = CliRunner()

    with runner.isolated_filesystem():
        result = runner.invoke(
            cli, ["alignments", "create-table", str(haplotagged_bam), "align_table.parquet", "--phased"]
        )
        assert result.exit_code == 0

        df = pd.read_parquet("align_table.parquet")
        print(df.haplotype.value_counts())


def test_reformat_bam(read_sorted_bam):
    runner = CliRunner()

    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["alignments", "reformat-bam", str(read_sorted_bam), "reformatted.sam"])
        assert result.exit_code == 0


def test_prepare_reads(read_fastq_file):
    runner = CliRunner()

    with runner.isolated_filesystem():
        cwd = Path(os.getcwd())
        result = runner.invoke(
            cli,
            [
                "reads",
                "prepare",
                str(read_fastq_file),
                "reads",
                "--min-qscore",
                "6",
                "--max-qscore",
                "12",
                "--min-read-length",
                "500",
                "--max-read-length",
                "5000",
            ],
        )
        assert result.exit_code == 0
        new_files = set([f.name for f in cwd.glob("*.*")])
        expected_files = set(
            [
                "reads.pass.fq.gz",
                "reads.fail.fq.gz",
                "reads.read_metadata.parquet",
                "reads.summary.csv",
                "reads.catalog.yaml",
            ]
        )

        df = pd.read_parquet(cwd / "reads.read_metadata.parquet", engine="pyarrow")
        # g = df.groupby("pass_filter")[['read_length', 'qscore']].agg(["min", "max", "count"])
        res = df["pass_filter"].value_counts()
        assert res[True] == 19
        assert new_files == expected_files


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

        result = runner.invoke(cli, ["refgenome", "fragments-to-hicref", "digest.fragments.parquet", "digest.hicRef"])
        assert result.exit_code == 0
        assert new_files == expected_files
