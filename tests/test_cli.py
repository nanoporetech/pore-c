import pandas as pd
from click.testing import CliRunner

from pore_c.cli import cli


def _run_command(opts):
    runner = CliRunner()
    common_opts = ["-vvv", "--dask-use-threads", "--dask-disable-dashboard", "--dask-scheduler-port", "0"]
    comd = common_opts + [str(_) for _ in opts]
    print("Running command: {}".format(" ".join(comd)))
    result = runner.invoke(cli, comd)
    return result


def test_parquet_to_csv(pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("utils")

    prefix = pore_c_table_pq.name.split(".")[0]

    output_csv = outdir / (prefix + ".pore_c.csv.gz")
    result = _run_command(["utils", "parquet-to-csv", pore_c_table_pq, output_csv])

    assert result.exit_code == 0


def test_haplotype_consensus(pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("haplotype_consensus")
    prefix = pore_c_table_pq.name.split(".")[0]

    output_pq = outdir / (prefix + ".pore_c.parquet")
    result = _run_command(
        ["alignments", "assign-consensus-haplotype", pore_c_table_pq, output_pq, "--threshold", "0.51"]
    )

    assert result.exit_code == 0
    pre_df = pd.read_parquet(pore_c_table_pq)[["align_idx", "haplotype"]].set_index("align_idx").sort_index()
    post_df = pd.read_parquet(output_pq)[["align_idx", "haplotype"]].set_index("align_idx").sort_index()
    # 6 alignments have their haplotypes changed as a result
    changes = (pre_df != post_df)["haplotype"].sum()
    assert changes == 7


def test_contacts_to_salsa_bed(contact_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(["contacts", "export", contact_table_pq, "salsa_bed", outdir / prefix])

    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = {prefix + ".salsa.bed"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_contacts_to_pairs(contact_table_pq, chromsizes, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        ["contacts", "export", contact_table_pq, "pairs", outdir / prefix, "--chromsizes", chromsizes]
    )

    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = {prefix + ".pairs"}
    assert result.exit_code == 0
    assert new_files == expected_files


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


def test_contacts_to_cool_by_haplotype(contact_table_pq, fragment_table_pq, chromsizes, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("to_cool_by_haplotype")
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
            "--by-haplotype",
        ]
    )
    assert result.exit_code == 0
    new_files = set([f.name for f in outdir.glob("*.*")])
    expected_files = {
        f"{prefix}.{ht_key}.cool" for ht_key in ("1_1", "1_2", "2_2", "nohap_1", "nohap_2", "nohap_nohap")
    }
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


def test_fragments_to_contacts(pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = pore_c_table_pq.name.split(".")[0]

    contact_table = str(outdir / (prefix + ".contacts.parquet"))

    result = _run_command(["alignments", "to-contacts", pore_c_table_pq, contact_table])
    assert result.exit_code == 0


def test_contact_summary(contact_table_pq, read_summary_csv, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    concatemer_table = str(outdir / (prefix + ".concatemer.parquet"))
    concatemer_summary = str(outdir / (prefix + ".concatemer.summary.csv"))

    result = _run_command(
        ["contacts", "summarize", contact_table_pq, read_summary_csv, concatemer_table, concatemer_summary]
    )
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


def test_create_table_with_haplotagged_aligns(coord_sorted_bam, haplotagged_aligns, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("create_table")
    align_table = str(outdir / "align_table.parquet")

    result = _run_command(
        [
            "alignments",
            "create-table",
            str(coord_sorted_bam),
            align_table,
            "--alignment-haplotypes",
            str(haplotagged_aligns),
        ]
    )
    assert result.exit_code == 0


def test_reformat_bam(read_sorted_bam, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("reformat_bam")
    result = _run_command(["alignments", "reformat-bam", str(read_sorted_bam), str(outdir / "reformatted.bam")])
    assert result.exit_code == 0


def test_prepare_reads(read_fastq_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("prepare_reads")

    result = _run_command(
        [
            "reads",
            "prepare",
            str(read_fastq_file),
            outdir / "reads",
            "--batch-size",
            "50",
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
            "reads.batch1.fq.gz",
            "reads.batch2.fq.gz",
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
