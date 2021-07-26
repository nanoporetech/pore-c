import pandas as pd
import pytest
from click.testing import CliRunner

from pore_c.cli import cli


@pytest.fixture(
    scope="session",
    params=[
        ["--dask-scheduler", "default"],
        pytest.param(["--dask-scheduler", "processes", "--dask-num-workers", "2"], marks=pytest.mark.skip),
        ["--dask-scheduler", "threads", "--dask-num-workers", "2"],
        pytest.param(
            [
                "--dask-scheduler",
                "local_cluster",
                "--dask-num-workers",
                "2",
                "--dask-threads-per-worker",
                "2",
                "--dask-disable-dashboard",
            ],
            marks=pytest.mark.skip,
        ),
    ],
    ids=["sched-default", "sched-processes", "sched-threads", "sched-local-cluster"],
)
def dask_settings(request):
    return request.param


def _run_command(opts, dask_settings=None):
    runner = CliRunner()
    common_opts = ["-vvv"]

    if dask_settings is None:
        dask_settings = ["--dask-scheduler", "default", "--dask-num-workers", "1"]
    comd = common_opts + dask_settings + [str(_) for _ in opts]
    print("Running command: {}".format(" ".join(comd)))
    result = runner.invoke(cli, comd)
    if result.exit_code != 0:
        print(result.output)
    return result


def test_parquet_to_csv(dask_settings, pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("utils")

    prefix = pore_c_table_pq.name.split(".")[0]

    output_csv = outdir / (prefix + ".pore_c.csv.gz")
    result = _run_command(["utils", "parquet-to-csv", pore_c_table_pq, output_csv], dask_settings=dask_settings)

    assert result.exit_code == 0


def test_haplotype_consensus(dask_settings, pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("haplotype_consensus")
    prefix = pore_c_table_pq.name.split(".")[0]

    output_pq = outdir / (prefix + ".pore_c.parquet")
    result = _run_command(
        ["alignments", "assign-consensus-haplotype", pore_c_table_pq, output_pq, "--threshold", "0.51"],
        dask_settings=dask_settings,
    )

    assert result.exit_code == 0
    pre_df = pd.read_parquet(pore_c_table_pq)[["align_idx", "haplotype"]].set_index("align_idx").sort_index()
    post_df = pd.read_parquet(output_pq)[["align_idx", "haplotype"]].set_index("align_idx").sort_index()
    # 6 alignments have their haplotypes changed as a result
    changes = (pre_df != post_df)["haplotype"].sum()
    assert changes == 7


def test_contacts_to_salsa_bed(dask_settings, contact_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        ["contacts", "export", contact_table_pq, "salsa_bed", outdir / prefix], dask_settings=dask_settings
    )

    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {prefix + ".salsa.bed"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_contacts_to_pairs(dask_settings, contact_table_pq, chromsizes, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        ["contacts", "export", contact_table_pq, "pairs", outdir / prefix, "--chromsizes", chromsizes],
        dask_settings=dask_settings,
    )

    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {prefix + ".pairs"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_contacts_to_paired_end_fastq(dask_settings, contact_table_pq, raw_refgenome_file, tmp_path_factory):
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
        ],
        dask_settings=dask_settings,
    )

    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {prefix + ".1.fastq", prefix + ".2.fastq"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_contacts_to_cool_by_haplotype(
    dask_settings, contact_table_pq, fragment_table_pq, chromsizes, tmp_path_factory
):
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
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0
    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {
        f"{prefix}.{ht_key}.cool" for ht_key in ("1_1", "1_2", "2_2", "nohap_1", "nohap_2", "nohap_nohap")
    }
    assert new_files == expected_files


def test_contacts_to_cool(dask_settings, contact_table_pq, fragment_table_pq, chromsizes, tmp_path_factory):
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
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_fragments_to_contacts(dask_settings, pore_c_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = pore_c_table_pq.name.split(".")[0]

    contact_table = str(outdir / (prefix + ".contacts.parquet"))

    result = _run_command(["alignments", "to-contacts", pore_c_table_pq, contact_table], dask_settings)
    assert result.exit_code == 0


def test_contact_summary(dask_settings, contact_table_pq, read_summary_csv, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    concatemer_table = str(outdir / (prefix + ".concatemer.parquet"))
    concatemer_summary = str(outdir / (prefix + ".concatemer.summary.csv"))

    result = _run_command(
        ["contacts", "summarize", contact_table_pq, read_summary_csv, concatemer_table, concatemer_summary],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_assign_fragments(dask_settings, align_table_pq, fragment_table_pq, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("assign_fragments")
    prefix = fragment_table_pq.name.split(".")[0]

    result = _run_command(
        [
            "alignments",
            "assign-fragments",
            str(align_table_pq),
            str(fragment_table_pq),
            str(outdir / (prefix + ".pore_c.parquet")),
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_create_table(dask_settings, haplotagged_bam, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("create_table")
    align_table = str(outdir / "align_table.parquet")

    result = _run_command(
        ["alignments", "create-table", str(haplotagged_bam), align_table], dask_settings=dask_settings
    )
    assert result.exit_code == 0


def test_create_table_with_haplotagged_aligns(dask_settings, coord_sorted_bam, haplotagged_aligns, tmp_path_factory):
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
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_reformat_bam(dask_settings, read_sorted_bam, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("reformat_bam")
    result = _run_command(
        ["alignments", "reformat-bam", str(read_sorted_bam), str(outdir / "reformatted.bam")],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_prepare_reads(dask_settings, read_fastq_file, tmp_path_factory):
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
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0
    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {
        "reads.batch1.fq.gz",
        "reads.batch2.fq.gz",
        "reads.fail.fq.gz",
        "reads.read_metadata.parquet",
        "reads.summary.csv",
        "reads.catalog.yaml",
    }
    df = pd.read_parquet(str(outdir / "reads.read_metadata.parquet"), engine="pyarrow")
    # g = df.groupby("pass_filter")[['read_length', 'qscore']].agg(["min", "max", "count"])
    res = df["pass_filter"].value_counts()
    assert res[True] == 19
    assert new_files == expected_files


def test_prepare_refgenome(dask_settings, raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("prepare_refgenome")
    result = _run_command(
        ["refgenome", "prepare", str(raw_refgenome_file), outdir / "refgenome"], dask_settings=dask_settings
    )
    assert result.exit_code == 0
    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {
        "refgenome.chromsizes",
        "refgenome.metadata.csv",
        "refgenome.fa.fai",
        "refgenome.catalog.yaml",
        "refgenome.fa",
    }
    assert new_files == expected_files


def test_virtual_digest(dask_settings, raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("virtual_digest")

    refgenome_prefix = raw_refgenome_file.name.split(".")[0]
    enzyme = "NlaIII"
    digest_id = f"{refgenome_prefix}_{enzyme}"

    result = _run_command(["refgenome", "prepare", str(raw_refgenome_file), outdir / refgenome_prefix])
    result = _run_command(
        ["refgenome", "virtual-digest", outdir / f"{refgenome_prefix}.fa", enzyme, outdir / digest_id], dask_settings
    )
    assert result.exit_code == 0
    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {
        f"{digest_id}.catalog.yaml",
        f"{digest_id}.digest_stats.csv",
        f"{digest_id}.fragments.parquet",
        f"{refgenome_prefix}.chromsizes",
        f"{refgenome_prefix}.metadata.csv",
        f"{refgenome_prefix}.fa.fai",
        f"{refgenome_prefix}.catalog.yaml",
        f"{refgenome_prefix}.fa",
    }
    assert new_files == expected_files

    result = _run_command(
        [
            "refgenome",
            "fragments-to-hicref",
            outdir / f"{digest_id}.fragments.parquet",
            outdir / f"{refgenome_prefix}.hicRef",
        ],
        dask_settings=dask_settings,
    )
    assert result.exit_code == 0


def test_contacts_to_merged_no_dups(dask_settings, contact_table_pq, raw_refgenome_file, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("contacts")
    prefix = contact_table_pq.name.split(".")[0]

    result = _run_command(
        [
            "contacts",
            "export",
            contact_table_pq,
            "merged_no_dups",
            outdir / prefix,
            "--reference-fasta",
            raw_refgenome_file,
        ],
        dask_settings=dask_settings,
    )

    new_files = {f.name for f in outdir.glob("*.*")}
    expected_files = {prefix + ".mnd.txt"}
    assert result.exit_code == 0
    assert new_files == expected_files


def test_extract_snv_links(
    dask_settings, coord_sorted_bam, phased_vcf, raw_refgenome_file_decompressed, tmp_path_factory
):
    outdir = tmp_path_factory.mktemp("extract_snv_links")
    prefix = coord_sorted_bam.name.split(".")[0]

    for x in range(2):
        result = _run_command(
            [
                "variants",
                "extract-snv-links",
                coord_sorted_bam,
                phased_vcf,
                raw_refgenome_file_decompressed,
                outdir / f"{prefix}_{x}.var_links.pq",
            ],
            dask_settings=dask_settings,
        )
        assert result.exit_code == 0

    new_files = {f.name for f in outdir.glob("*.*")}

    dfs = [pd.read_parquet(outdir / f) for f in new_files]

    agg_counts = pd.concat([_[["cis_count", "trans_count"]] for _ in dfs], axis=0).sum().to_dict()
    expected_files = {f"{prefix}_{x}.var_links.pq" for x in range(2)}
    result = _run_command(
        [
            "variants",
            "aggregate-links",
            outdir / f"{prefix}_0.var_links.pq",
            outdir / f"{prefix}_1.var_links.pq",
            outdir / f"{prefix}_agg.var_links.pq",
        ],
        dask_settings=dask_settings,
    )

    merged_df = pd.read_parquet(outdir / f"{prefix}_agg.var_links.pq")
    assert agg_counts == merged_df[["cis_count", "trans_count"]].sum().to_dict()
    assert result.exit_code == 0

    assert new_files == expected_files
