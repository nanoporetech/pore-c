import pandas as pd
from click.testing import CliRunner

from pore_c.cli import cli


def _run_command(opts):
    runner = CliRunner()
    common_opts = ["-vvv", "--dask-use-threads", "--dask-disable-dashboard", "--dask-scheduler-port", "0"]
    comd = common_opts + [str(_) for _ in opts]
    # print("Running command: {}".format(" ".join(comd)))
    result = runner.invoke(cli, comd)
    return result


def test_contacts_to_paired_end_fastq(contact_table_pq, raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        [
            "contacts",
            "export",
            contact_table_pq,
            "paired_end_fastq",
            outdir / prefix,
            "--reference-fasta",
            raw_refgenome_file,
        ]
    )

    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = {prefix + ".1.fastq", prefix + ".2.fastq"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_contacts_to_cool(contact_table_pq, fragment_table_pq, chromsizes, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        [
            "contacts",
            "export",
            contact_table_pq,
            "cooler",
            outdir / prefix,
            "--fragment-table",
            fragment_table_pq,
            "--chromsizes",
            chromsizes,
        ]
    )
    assert result.exit_code == 0


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

    refgenome_prefix = raw_refgenome_file.name.split(".")[0]
    enzyme = "NlaIII"
    digest_id = f"{refgenome_prefix}_{enzyme}"

    result = _run_command(["refgenome", "prepare", str(raw_refgenome_file), outdir / refgenome_prefix])
    result = _run_command(
        ["refgenome", "virtual-digest", outdir / f"{refgenome_prefix}.fa", enzyme, outdir / digest_id]
    )
    assert result.exit_code == 0
    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = set(
        [
            f"{digest_id}.catalog.yaml",
            f"{digest_id}.digest_stats.csv",
            f"{digest_id}.fragments.parquet",
            f"{refgenome_prefix}.chromsizes",
            f"{refgenome_prefix}.metadata.csv",
            f"{refgenome_prefix}.fa.fai",
            f"{refgenome_prefix}.catalog.yaml",
            f"{refgenome_prefix}.fa",
        ]
    )
    assert new_files == expected_files

    result = _run_command(
        [
            "refgenome",
            "fragments-to-hicref",
            outdir / f"{digest_id}.fragments.parquet",
            outdir / f"{refgenome_prefix}.hicRef",
        ]
    )
    assert result.exit_code == 0
