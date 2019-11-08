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
There are some sample datasets (fastq, alignment parquets, .pairs, .cool files) available for HindIII-digested HG002 (31Gb) [here](https://ont-applications-aws-self-service-us-east-1.s3.amazonaws.com/0aface03-968e-42b8-a0c0-c81068f11206/porec_sync/20191103.preprint_HG002.tar.gz?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJ7PQ6GQTKDJTLDRQ%2F20191108%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20191108T175453Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=b66e41ea9c704c1bebb079b46b43085be4a8a3552a4697a3db43e888b5da6efc) and for NlaIII-digested GM12878 (23Gb) [here](https://ont-applications-aws-self-service-us-east-1.s3.amazonaws.com/0aface03-968e-42b8-a0c0-c81068f11206/porec_sync/20191103.preprint_NA12878.tar.gz?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJ7PQ6GQTKDJTLDRQ%2F20191108%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20191108T175802Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=281f4c581855960f998da3fbe4e3ba670b509fc24f792a64f2cc079def5c1ea1).

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
