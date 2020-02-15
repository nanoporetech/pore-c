import pandas as pd
from click.testing import CliRunner

from pore_c.cli import cli


def _run_command(opts):
    runner = CliRunner()
    common_opts = ["-vvv", "--dask-use-threads", "--dask-disable-dashboard", "--dask-scheduler-port", "0"]

    result = runner.invoke(cli, common_opts + [str(_) for _ in opts])
    return result


def test_fragments_to_contacts(align_table_pq, fragment_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = fragment_table_pq.name.split(".")[0]

    fragment_table = str(outdir / (prefix + ".pore_c.parquet"))
    contact_table = str(outdir / (prefix + ".contacts.parquet"))
    concatemer_table = str(outdir / (prefix + ".concatemers.parquet"))

    result = _run_command(
        ["alignments", "assign-fragments", str(align_table_pq), str(fragment_table_pq), fragment_table]
    )
    assert result.exit_code == 0
    result = _run_command(["alignments", "to-contacts", fragment_table, contact_table, concatemer_table])
    assert result.exit_code == 0


def test_assign_fragments(align_table_pq, fragment_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("assign_fragments")
    prefix = fragment_table_pq.name.split(".")[0]

    result = _run_command(
        [
            "alignments",
            "assign-fragments",
            str(align_table_pq),
            str(fragment_table_pq),
            str(outdir / (prefix + ".pore_c.parquet")),
        ]
    )
    assert result.exit_code == 0


def test_create_table(haplotagged_bam, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("create_table")
    align_table = str(outdir / "align_table.parquet")

    result = _run_command(["alignments", "create-table", str(haplotagged_bam), align_table])
    assert result.exit_code == 0
    df = pd.read_parquet(align_table)
    print(df.haplotype.value_counts())


def test_reformat_bam(read_sorted_bam, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("reformat_bam")

    result = _run_command(["alignments", "reformat-bam", str(read_sorted_bam), str(outdir / "reformatted.sam")])
    assert result.exit_code == 0


def test_prepare_reads(read_fastq_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("prepare_reads")

    result = _run_command(
        [
            "reads",
            "prepare",
            str(read_fastq_file),
            outdir / "reads",
            "--min-qscore",
            "6",
            "--max-qscore",
            "12",
            "--min-read-length",
            "500",
            "--max-read-length",
            "5000",
        ]
    )
    assert result.exit_code == 0
    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = set(
        [
            "reads.pass.fq.gz",
            "reads.fail.fq.gz",
            "reads.read_metadata.parquet",
            "reads.summary.csv",
            "reads.catalog.yaml",
        ]
    )
    df = pd.read_parquet(str(outdir / "reads.read_metadata.parquet"), engine="pyarrow")
    # g = df.groupby("pass_filter")[['read_length', 'qscore']].agg(["min", "max", "count"])
    res = df["pass_filter"].value_counts()
    assert res[True] == 19
    assert new_files == expected_files


def test_prepare_refgenome(raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("prepare_refgenome")
    result = _run_command(["refgenome", "prepare", str(raw_refgenome_file), outdir / "refgenome"])
    assert result.exit_code == 0
    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = set(
        ["refgenome.chromsizes", "refgenome.metadata.csv", "refgenome.fa.fai", "refgenome.catalog.yaml", "refgenome.fa"]
    )
    assert new_files == expected_files


def test_virtual_digest(raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("virtual_digest")

    result = _run_command(["refgenome", "prepare", str(raw_refgenome_file), outdir / "refgenome"])
    result = _run_command(["refgenome", "virtual-digest", outdir / "refgenome.fa", "NlaIII", outdir / "digest"])
    assert result.exit_code == 0
    new_files = set([f.name for f in outdir.glob("*.*")])
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

    result = _run_command(
        ["refgenome", "fragments-to-hicref", outdir / "digest.fragments.parquet", outdir / "digest.hicRef"]
    )
    assert result.exit_code == 0
