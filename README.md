# concatemer_finding
Flags sequences that are entirely composed of tandem repeats

## Install:

Navigate to the repositorie's root dir:

    cd concatemer_finding

Create a conda envirnoment:

    mamba env create -p ./conda.env -f conda.yaml

Activate the environment

    conda activate ./conda.env

Just run snakemake to teste

    snakemake

## Run

Change the default configuration for snakemake. This assumes you run the repo
dir. See the snakemake docs for running elsewhere.

    snakemake --config output_dir='/path/to/output_dir' \
                       reads_template='/path/to/{sample}.reads.fasta'

The above assumes reads from each sample are in a single fasta file per sample.
