# Command Line Interface

```{contents} Table of Contents
---
depth: 3
---
```

## refgenome

### prepare

Pre-process a reference genome

#### Usage

```bash
pore_c refgenome prepare [OPTIONS] REFERENCE_FASTA OUTPUT_PREFIX

Pre-process a reference genome for use by pore-C tools.

Prepare a reference genome or draft assembly use by pore-c tools.
This tool creates the following files:

    <output_prefix>.fa - A decompressed fasta file
    <output_prefix>.chromsizes - A tab-separated list of chromosome lengths
    <output_prefix>.metadata.csv - A csv with metadata about each of
    <output_prefix>.catalog.yaml - An intake catalog of these files
```

#### Parameters:

-   _reference_fasta_ [required]
-   _output_prefix_ [required]
-   _--genome-id TEXT_: An ID for this genome assembly

---

### virtual-digest

Virtual digest of a reference genome.

#### Usage

```bash
pore_c refgenome virtual-digest [OPTIONS] FASTA CUT_ON OUTPUT_PREFIX

Carry out a virtual digestion of the genome/assembly FASTA.

The DIGEST_TYPE sets what type of value you can use for CUT_ON.


    ---------------------------------------------------------------------
    | digest_type |  cut_on           | notes                           |
    ---------------------------------------------------------------------
    | enzyme      | NlaIII            | Enzyme name is case sensitive   |
    | enzyme      | HindIII           |                                 |
    | regex       | (GAATTC|GCGGCCGC) | Two site EcoRI and NotI         |
    | regex       | RAATY             | Degenerate site ApoI            |
    | bin         | 50000             | Create equal-width bins of 50k  |
    | bin         | 50k               | Create equal-width bins of 50k  |
    =====================================================================

This tool will create the following output files:


    <output_prefix>.fragments.parquet

      A table containing the coordinates and some metadata on each fragment

    <output_prefix>.digest_stats.csv

      Per chromosome/contig summary statistics

    <output_prefix>.catalog.yaml

      An intake catalog of the digest files
```

#### Parameters:

-   _fasta_ [required]
-   _cut_on_ [required]
-   _output_prefix_ [required]
-   _--digest-type [enzyme|regex|bin]_: The type of digest you want to do
-   _-n, --n_workers INTEGER_: The number of dask_workers to use [default: 1]

---

### fragments-to-hicref

Create a hicRef file for a virtual digest.

#### Usage

```bash
pore_c refgenome fragments-to-hicref [OPTIONS] FRAGMENTS_PARQUET HICREF

Convert a  .fragments.parquet file to hicRef format
```

#### Parameters:

-   _fragments_parquet_ [required]
-   _hicref_ [required]

---

## reads

### prepare

Create a catalog file for a set of reads

#### Usage

```bash
pore_c reads prepare [OPTIONS] FASTQ OUTPUT_PREFIX

Preprocess a set of reads for use with pore_c tools.

This tool creates the following files:


    <output_prefix>.batch[batch_idx].fq.gz - Fastq with all the reads that pass the qscore and
      length filters. Fastqs are split so there are at most --batch-size reads per fastq
    <output_prefix>.fail.fq.gz - Reads that fail the filters
    <output_prefix>.read_metadata.parquet - Length and qscore metadata for each read and
      whether they pass the filter
    <output_prefix>.summary.csv - Summary stats for all/pass/fail reads
    <output_prefix>.catalog.yaml - An intake catalog
```

#### Parameters:

-   _fastq_ [required]
-   _output_prefix_ [required]
-   _--batch-size INTEGER_: The reads will be split into batches of this size for downstream processing [default: 10000]
-   _--min-read-length INTEGER_: The minimum length read to run through pore_c [default: 1]
-   _--max-read-length INTEGER_: The maximum length read to run through pore_c. Note that bwa mem can crash on very long reads [default: 150000]
-   _--min-qscore INTEGER_: The minimum read qscore [default: 0]
-   _--max-qscore INTEGER_: The maximum read qscore [default: 266]
-   _--user-metadata TEXT_: Additional user metadata to associate with this run

---

## alignments

### reformat-bam

Reformat a BAM file to have a unique read name per alignment

#### Usage

```bash
pore_c alignments reformat-bam [OPTIONS] INPUT_SAM OUTPUT_SAM

Reformat query_name in INPUT_SAM  and write to OUTPUT_SAM

This tool reformats an alignment file so that it works with downstream
steps in the Pore-C pipeline. For both files you can supply '-' if you want
to read/write from/to stdin/stdout. The 'query_name' field of the alignment
file will be reformatted so that each alignment in the SAM file has a
unique query name:


    <read_id> -> <read_id>:<read_idx>:<align_idx>

Where 'read_idx' is a unique integer id for each read within the file and
'align_idx' is a unique integer id for each alignment within the file. The
tool also adds a 'BX' tag consisting of the 'read_id' to each record.
```

#### Parameters:

-   _input_sam_ [required]
-   _output_sam_ [required]
-   _--input-is-bam_: If piping a BAM from stdin (rather than sam) [default: False]
-   _--output-is-bam_: If piping a BAM to stdout (rather than sam) [default: False]
-   _--set-bx-flag_: Set the BX tag to the read name [default: False]

---

### create-table

Parse a namesortd bam to pore-C alignment format

#### Usage

```bash
pore_c alignments create-table [OPTIONS] INPUT_BAM OUTPUT_TABLE

Convert a BAM file to a tabular format sorted by read for downstream analysis
```

#### Parameters:

-   _input_bam_ [required]
-   _output_table_ [required]
-   _--alignment-haplotypes PATH_: The alignment to haplotype mapping from whatshap

---

### assign-fragments

Parse a namesortd bam to pore-C alignment format

#### Usage

```bash
pore_c alignments assign-fragments [OPTIONS] ALIGN_TABLE FRAGMENTS_TABLE PORE_C_TABLE

For each alignment in ALIGN_TABLE either filter out or assign a fragment from FRAGMENT_TABLE
```

#### Parameters:

-   _align_table_ [required]
-   _fragments_table_ [required]
-   _pore_c_table_ [required]
-   _--mapping_quality_cutoff INTEGER_: Minimum mapping quality for an alignment to be considered [default: 1]
-   _--min_overlap_length INTEGER_: Minimum overlap in base pairs between an alignment and restriction fragment [default: 10]
-   _--containment_cutoff FLOAT_: Minimum percentage of a fragment included in an overlap for that fragment to be considered 'contained' within an alignment [default: 99.0]

---

### filter-bam

Filter bam using pore_c table

#### Usage

```bash
pore_c alignments filter-bam [OPTIONS] INPUT_BAM PORE_C_TABLE OUTPUT_BAM


```

#### Parameters:

-   _input_bam_ [required]
-   _pore_c_table_ [required]
-   _output_bam_ [required]
-   _--clean-read-name_: Strip out the extra information placed in the BAM by reformat_bam

---

### assign-consensus-haplotype

Parse a namesortd bam to pore-C alignment format

#### Usage

```bash
pore_c alignments assign-consensus-haplotype [OPTIONS] PORE_C_TABLE OUTPUT_PORE_C_TABLE

Calculate a per-read consensus haplotype for each phase_set in ALIGN_TABLE and write the results back
to OUTPUT_ALIGN_TABLE
```

#### Parameters:

-   _pore_c_table_ [required]
-   _output_pore_c_table_ [required]
-   _--threshold FLOAT_: major:minor haplotype fraction must be greater than this value to assign a consensus [default: 0.8]

---

### to-contacts

Parses the alignment table and converts to pairwise contacts

#### Usage

```bash
pore_c alignments to-contacts [OPTIONS] PORE_C_TABLE CONTACT_TABLE

Covert the alignment table to a pairwise contact table and associated concatemer table
```

#### Parameters:

-   _pore_c_table_ [required]
-   _contact_table_ [required]

---

## contacts

### merge

Summarise a contact table

#### Usage

```bash
pore_c contacts merge [OPTIONS] [SRC_CONTACT_TABLES]... DEST_CONTACT_TABLE


```

#### Parameters:

-   _src_contact_tables_ [required]
-   _dest_contact_table_ [required]
-   _--fofn_: If this flag is set then the SRC_CONTACT_TABLES is a file of filenames corresponding to the contact tables you want to merge. This is workaround for when the command line gets too long.

---

### downsample

Downsample a contact table

#### Usage

```bash
pore_c contacts downsample [OPTIONS] SRC_CONTACT_TABLE DEST_CONTACT_TABLE_PREFIX
        [DOWNSAMPLE_INCREMENTS]...


```

#### Parameters:

-   _src_contact_table_ [required]
-   _dest_contact_table_prefix_ [required]
-   _downsample_increments_ [required]
-   _--downsample-unit [Gb|Mb|Kb]_:
-   _--random-seed INTEGER_:
-   _--tol FLOAT_: Check if the difference between the sampled amout and the target amount is greater than this proportion
-   _--warn_: If the a sample fails the --tol check print a warning rather than exiting
-   _--max-attempts INTEGER_: The number of times to try and find a set of subsamples all within --tol

---

### summarize

Summarise a contact table

#### Usage

```bash
pore_c contacts summarize [OPTIONS] CONTACT_TABLE READ_SUMMARY_TABLE CONCATEMER_TABLE
        CONCATEMER_SUMMARY_CSV


```

#### Parameters:

-   _contact_table_ [required]
-   _read_summary_table_ [required]
-   _concatemer_table_ [required]
-   _concatemer_summary_csv_ [required]
-   _--user-metadata TEXT_: Add additional user metadata to the summary table, must be a dictionary in json format

---

### export

Export contacts to various formats

#### Usage

```bash
pore_c contacts export [OPTIONS] CONTACT_TABLE
        [cooler|salsa_bed|paired_end_fastq|pairs|merged_no_dups] OUTPUT_PREFIX

Export contacts to the following formats:

 - cooler: a sparse representation of a contact matrix
 - paired_end_fastq: for each contact create a pseudo pair-end read using the reference genome sequence
```

#### Parameters:

-   _contact_table_ [required]
-   _format_ [required]
-   _output_prefix_ [required]
-   _--min-mapping-quality INTEGER_: Both alignments have mapping qualities greater than this [default: 0]
-   _--min-align-base-qscore INTEGER_: Both alignments have mean base qualities greater than this [default: 0]
-   _--cooler-resolution INTEGER_: The bin width of the resulting matrix [default: 1000]
-   _--fragment-table TEXT_: The fragment table for the corresponding virtual digest(required if export format is in cooler)
-   _--by-haplotype_: Create a cooler for each pair of haplotypes (eg 1-1, 1-2, 2-2,...). Only valid with 'cooler'
-   _--chromsizes TEXT_: The chromsizes file for the corresponding genome(required if export format is in cooler,pairs)
-   _--reference-fasta TEXT_: The reference genome used to generate the contact table(required if export format is in paired_end_fastq,merged_no_dups)

---

## utils

### parquet-to-csv

#### Usage

```bash
pore_c utils parquet-to-csv [OPTIONS] INPUT_PARQUET OUTPUT_CSV

Convert a parquet file to CSV
```

#### Parameters:

-   _input_parquet_ [required]
-   _output_csv_ [required]
