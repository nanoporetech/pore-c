![.](.docs/_static/ONT_logo.png "Oxford Nanopore Technologies")

# Pore-C tools

<!-- start pore-c-overview -->

This package is designed to analyse the data from multi-contact pore-C reads. It is similar to the
[pairtools](https://github.com/mirnylab/pairtools) package in scope, however it is specifically designed to
handle multi-contact reads (aka c-walks). It carries out the following
operations:

-   pre-processing a reference genome to generate auxiliary files used in downstream analyses
-   creating virtual digests of the reference genome
-   processing BAM files to filter spurious alignments, detect ligation junctions and assign fragments
-   converting the resulting contacts to a [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) and a COO-formatted matrix compatible with [Cooler](https://github.com/mirnylab/cooler) for downstream processing.

The recommended way to run Pore-C tools is through the associated [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake) pipeline. It wraps the `pore-C` commands and also handles some of the analysis steps outside the scope of `pore-C tools` such as read alignment and conversion of output files to `.cool` format. If you want to run
the tools outside of this pipeline please see the section below on

<!-- end pore-c-overview -->

### Sample data

There are some sample datasets (fastq, alignment parquets, .pairs, .cool files) available for HindIII-digested HG002 (31Gb) [here](https://ont-datasets-us-east-1-public.s3.amazonaws.com/20191103.preprint_HG002.tar.gz) and for NlaIII-digested GM12878 (23Gb) [here](https://ont-datasets-us-east-1-public.s3.amazonaws.com/20191103.preprint_NA12878.tar.gz).

## Development

The recommended way to use pore-C tools is through the associated [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake) pipeline, however if you want to run tools outside of this context you'll need to set
up a development environment. We use the [tox](https://tox.readthedocs.io/en/latest/index.html) automation tool with the [tox-conda](https://github.com/tox-dev/tox-conda) extension
to create conda environments for development. Several environments are defined in the `tox.ini` file:

-   **py37:** The default environment. Running `tox -e py37` will run all of the `pytest` tests. By adding `--` to the command you you can pass additional arguments to pytest e.g. `tox -e py37 -- -k refgenome`.
-   **dev:** A development environment that runs the `pore_c` command. You can use this to try our `pore_c` commands e.g. `tox -e dev -- refgenome --help`.
-   **snakemake:** An environment that runs an associated snakemake pipeline in development mode. The command assumes that the `pore-c-snakemake` repository is located in the same directory as `pore-c` though this can be changed with enviroment variables - see the `tox.ini` file for details.

For each of these environments a conda environment is created at `.tox/<env_id>`. You can activate these for interactive work with `conda activate .tox/<env_id>`.

#### Developing the Snakemake pipeline

If you're working on the [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake) pipeline:

    git clone git@github.com:nanoporetech/pore-c.git
    cd pore-c
    tox -e snakemake --notest

    cd ../
    git clone git@github.com:nanoporetech/Pore-C-Snakemake.git
    cd pore-c-snakemake
    conda activate ../pore-c/.tox/snakemake
    snakemake -j 8 test --use-conda --config pore_c_version=dev

### Citation

The biorxiv pre-print describing Pore-C can be found here:

[Nanopore sequencing of DNA concatemers reveals higher-order features of chromatin structure](https://www.biorxiv.org/content/10.1101/833590v1)
