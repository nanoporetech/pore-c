import logging
from collections import defaultdict
from contextlib import ExitStack
from itertools import combinations
from pathlib import Path

import pandas as pd
import pysam
from whatshap.cli import CommandLineError, PhasedInputReader
from whatshap.cli.haplotag import get_variant_information, load_chromosome_variants
from whatshap.core import NumericSampleIds
from whatshap.vcf import VcfError, VcfReader

from pore_c.model import VariantPairRecord

logger = logging.getLogger(__name__)


def extract_snv_links(  # noqa: C901
    alignment_file: Path, variant_file: Path, reference: Path, mapping_quality_cutoff=10, variant_quality_cutoff=10
):

    with ExitStack() as stack:
        try:
            vcf_reader = stack.enter_context(VcfReader(variant_file, indels=True, phases=True))
        except OSError as err:
            raise CommandLineError(f"Error while loading variant file {variant_file}: {err}")

        use_vcf_samples = vcf_reader.samples
        try:
            bam_reader = stack.enter_context(pysam.AlignmentFile(alignment_file, "rb", require_index=True))
        except OSError as err:
            raise CommandLineError(f"Error while loading alignment file {alignment_file}: {err}")
        # This checks also sample compatibility with VCF
        shared_samples = use_vcf_samples

        # Check if user has specified a subset of regions per chromosome
        user_regions = {r: [(0, None)] for r in bam_reader.references}

        ignore_read_groups = True
        phased_input_reader = stack.enter_context(
            PhasedInputReader([alignment_file], reference, NumericSampleIds(), ignore_read_groups, indels=False)
        )

        data = defaultdict(set)
        for chrom, regions in user_regions.items():
            logger.debug(f"Processing chromosome {chrom}")

            # If there are no alignments for this chromosome, skip it. This allows to have
            # extra chromosomes in the BAM compared to the VCF as long as they are not actually
            # used.
            has_any_alignments = False
            for _ in bam_reader.fetch(contig=chrom):
                has_any_alignments = True
                break
            if not has_any_alignments:
                continue
            try:
                variant_table = load_chromosome_variants(vcf_reader, chrom, regions)
            except VcfError as e:
                raise CommandLineError(str(e))
            if variant_table is not None:
                logger.debug("Preparing haplotype information")
                for sample in shared_samples:
                    _, variants = get_variant_information(variant_table, sample)
                    read_set, _ = phased_input_reader.read(variant_table.chromosome, variants, sample, regions=regions)
                    for read in read_set:
                        if read.mapqs[0] < mapping_quality_cutoff:
                            continue
                        read_name = read.name.split(":")[0]
                        for var in read:
                            if var.quality < variant_quality_cutoff:
                                continue
                            data[read_name].add((chrom, var.position, var.allele))

        res = []
        for read_name, variants in data.items():
            for (var_1, var_2) in combinations(variants, 2):
                is_cis = var_1[-1] == var_2[-1]
                res.append((var_1[0], var_1[1], var_2[0], var_2[1], is_cis))

        df = pd.DataFrame(res, columns=["var1_chrom", "var1_position", "var2_chrom", "var2_position", "cis"])
        chrom_dtype = pd.CategoricalDtype(bam_reader.references, ordered=True)

        counts = (
            df.groupby(["var1_chrom", "var1_position", "var2_chrom", "var2_position", "cis"])
            .size()
            .unstack(fill_value=0)
            .rename(columns={True: "cis_count", False: "trans_count"})
            .reset_index()
            .astype(VariantPairRecord.pandas_dtype(overrides={"var1_chrom": chrom_dtype, "var2_chrom": chrom_dtype}))
            .sort_values(["var1_chrom", "var2_chrom", "var1_position", "var2_position"])
        )

        return counts
