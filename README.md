# Pore-C tools

This package is designed to analyse the data from multi-contact pore-C reads. It is similar to the
[pairtools](https://github.com/mirnylab/pairtools) package in scope, however it is specifically designed to
handle multi-contact reads (aka c-walks). It carries out the following
operations:

- pre-processing a reference genome to generate auxiliary files used in downstream analyses
- creating virtual digests of the reference genome
- processing read-sorted BAM files to filter spurious alignments, detect ligation junctions and assign fragments
- converting the resulting contacts to a [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) and a COO-formatted matrix compatible with [Cooler](https://github.com/mirnylab/cooler) for downstream processing.

There is an associated snakemake pipeline available at **TBD** that wraps the `pore-C` commands and also handles some of the analysis steps outside the scope of `pore-C tools` such as read alignment and conversion of output files to `.cool` format. This is the recommended way to run `pore-C tools`.

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
