# Pore-C tools

This package is designed to analyse the data from multi-contact pore-C reads. It is similar to the
[pairtools](https://github.com/mirnylab/pairtools) package in scope, however it is specifically designed to
handle multi-contact reads (aka c-walks). It carries out the following
operations:

- pre-processing a reference genome to generate auxiliary files used in downstream analyses
- creating virtual digests of the reference genome
- processing read-sorted BAM files to filter spurious alignments, detect ligation junctions and assign fragments
- converting the resulting contacts to a [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) and a COO-formatted matrix compatible with [Cooler](https://github.com/mirnylab/cooler) for downstream processing.

There is an associated [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake) that wraps the `pore-C` commands and also handles some of the analysis steps outside the scope of `pore-C tools` such as read alignment and conversion of output files to `.cool` format. This is the recommended way to run `pore-C tools`.

### Sample data
There are some sample datasets (fastq, alignment parquets, .pairs, .cool files) available for HindIII-digested HG002 (31Gb) [here](https://ont-datasets-us-east-1-public.s3.amazonaws.com/20191103.preprint_HG002.tar.gz) and for NlaIII-digested GM12878 (23Gb) [here](https://ont-datasets-us-east-1-public.s3.amazonaws.com/20191103.preprint_NA12878.tar.gz).

## Getting Started

### Creating a development environment

```
conda env create
source activate poreC
pip install -e .
pore_c --help
```

### Running tests

```
pytest tests
```

### Citation

The biorxiv pre-print describing Pore-C can be found here: 

[Nanopore sequencing of DNA concatemers reveals higher-order features of chromatin structure](https://www.biorxiv.org/content/10.1101/833590v1)
