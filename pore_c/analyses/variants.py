import logging
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple

import pandas as pd
import pysam
from whatshap.cli import CommandLineError, PhasedInputReader
from whatshap.cli.haplotag import get_variant_information, load_chromosome_variants
from whatshap.core import NumericSampleIds
from whatshap.vcf import VcfError, VcfReader

from pore_c.config import PQ_ENGINE
from pore_c.io import TableWriter
from pore_c.model import VariantPhaseRecord, VariantPhaseRecordDf

logger = logging.getLogger(__name__)


@dataclass
class QualCutoffs:
    mapping: int
    base: int


def extract_variant_info(  # noqa: C901
    output_parquet: Path,
    alignment_file: Path,
    variant_file: Path,
    reference: Path,
    mapping_quality_cutoff: int = 5,
    base_quality_cutoff: int = 5,
    pore_c_table: Optional[Path] = None,
) -> int:

    cutoffs = QualCutoffs(mapping_quality_cutoff, base_quality_cutoff)
    if pore_c_table:
        aligns = pd.read_parquet(pore_c_table, engine=PQ_ENGINE, columns=["align_idx", "pass_filter"]).set_index(
            ["align_idx"]
        )
        aligns = aligns[aligns["pass_filter"]].index
        if len(aligns) == 0:
            logger.warning("No aligments passed filter")
    else:
        aligns = None

    num_records = 0
    with ExitStack() as stack:
        try:
            vcf_reader = stack.enter_context(VcfReader(variant_file, indels=True, phases=True))
        except OSError as err:
            raise CommandLineError(f"Error while loading variant file {variant_file}: {err}")

        use_vcf_samples = vcf_reader.samples
        if use_vcf_samples is None:
            vcf_sample_id = None
        elif len(use_vcf_samples) > 1:
            raise NotImplementedError(f"Phased VCF must be single-sample: {variant_file} {use_vcf_samples}")
        else:
            vcf_sample_id = use_vcf_samples[0]
        try:
            bam_reader = stack.enter_context(pysam.AlignmentFile(alignment_file, "rb", require_index=True))
        except OSError as err:
            raise CommandLineError(f"Error while loading alignment file {alignment_file}: {err}")

        ignore_read_groups = True
        phased_input_reader = stack.enter_context(
            PhasedInputReader([alignment_file], reference, NumericSampleIds(), ignore_read_groups, indels=False)
        )

        writer = TableWriter(output_parquet)

        for table in iter_chromosomes(
            vcf_reader, bam_reader, phased_input_reader, vcf_sample_id, cutoffs, aligns=aligns
        ):
            # dump in batches by chromosome
            logger.debug(f"Writing {len(table)} records to {output_parquet}")
            writer(table)

        num_records = writer.row_counter
        writer.close()

    return num_records


#
#        res = []
#        for read_name, variants in data.items():
#            for (var_1, var_2) in combinations(variants, 2):
#                is_cis = var_1[-1] == var_2[-1]
#                res.append((var_1[0], var_1[1], var_2[0], var_2[1], is_cis))
#
#        df = pd.DataFrame(res, columns=["var1_chrom", "var1_position", "var2_chrom", "var2_position", "cis"])
#        chrom_dtype = pd.CategoricalDtype(bam_reader.references, ordered=True)
#
#        counts = (
#            df.groupby(["var1_chrom", "var1_position", "var2_chrom", "var2_position", "cis"])
#            .size()
#            .unstack(fill_value=0)
#            .rename(columns={True: "cis_count", False: "trans_count"})
#            .reset_index()
#        )
#
#        try:
#            counts = (counts
#            .astype(VariantPairRecord.pandas_dtype(overrides={"var1_chrom": chrom_dtype, "var2_chrom": chrom_dtype}))
#            .sort_values(["var1_chrom", "var2_chrom", "var1_position", "var2_position"])
#            )
#        except:
#            print(counts.head())
#            print(counts.dtypes)
#            raise
#
#        return counts


def iter_chromosomes(
    vcf_reader: VcfReader,
    bam_reader: pysam.AlignmentFile,
    phased_input_reader: PhasedInputReader,
    vcf_sample_id: Optional[str],
    cutoffs: QualCutoffs,
    aligns: Optional[Iterable] = None,
) -> Iterator[VariantPhaseRecordDf]:
    vcf_contigs = set(vcf_reader._vcf_reader.header.contigs)
    bam_contigs = set(bam_reader.references)
    has_alignments = set()
    for chrom in bam_contigs:
        for _ in bam_reader.fetch(contig=chrom):
            has_alignments.add(chrom)

    final_chroms = vcf_contigs & bam_contigs & has_alignments

    if final_chroms != vcf_contigs:
        logger.info(
            f"The following chromosomes in the VCF have no corresponding alignments: {vcf_contigs - final_chroms}"
        )

    if final_chroms != bam_contigs:
        logger.info(
            f"The following chromosomes in the BAM either have no corresponding alignments or no variants : {bam_contigs - final_chroms}"
        )

    chrom_dtype = pd.CategoricalDtype(bam_reader.references, ordered=True)

    for chrom in final_chroms:
        table = get_read_variants_by_chrom(
            chrom, vcf_reader, phased_input_reader, vcf_sample_id, cutoffs, aligns=aligns, chrom_dtype=chrom_dtype
        )
        if table is None:
            logger.warning(f"No variants found on {chrom}")
        else:
            yield table


def get_read_variants_by_chrom(  # noqa: C901
    chrom: str,
    vcf_reader: VcfReader,
    phased_input_reader: PhasedInputReader,
    vcf_sample_id: Optional[str],
    cutoffs: QualCutoffs,
    aligns: Optional[Iterable] = None,
    regions: Optional[List[Tuple[int, int]]] = None,
    chrom_dtype: Optional[pd.CategoricalDtype] = None,
) -> VariantPhaseRecordDf:

    if regions is None:
        regions = [(0, None)]
    try:
        variant_table = load_chromosome_variants(vcf_reader, chrom, regions)
    except VcfError as e:
        raise CommandLineError(str(e))
    if variant_table is None:
        logger.debug(f"No variants found for {chrom}")
        return None
    variantpos_to_phaseinfo, variants = get_variant_information(variant_table, vcf_sample_id)
    read_set, _ = phased_input_reader.read(variant_table.chromosome, variants, vcf_sample_id, regions=regions)

    results = []
    for read in read_set:
        try:
            read_name, _, align_idx = read.name.split(":")
            align_idx = int(align_idx)
        except ValueError:
            align_idx = None
            read_name = read.name
        if read.mapqs[0] < cutoffs.mapping:
            continue
        if aligns is not None:
            if align_idx is not None and align_idx not in aligns:
                continue
        read_name = read.name.split(":")[0]
        for var in read:
            base_qual = var.quality
            if base_qual < cutoffs.base:
                continue
            phaseset, haplotype = variantpos_to_phaseinfo[var.position]
            results.append(
                (
                    read_name,
                    align_idx,
                    chrom,
                    var.position,
                    var.allele,
                    read.mapqs[0],
                    base_qual,
                    phaseset,
                    haplotype,
                )
            )

    if len(results) == 0:
        return None

    dtype = VariantPhaseRecord.pandas_dtype(overrides={"chrom": chrom_dtype})
    names = dtype.keys()
    results = pd.DataFrame(results, columns=names).astype(dtype).sort_values(["read_name", "chrom", "position"])
    return results
