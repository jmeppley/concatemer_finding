# concatemer_finding
Flags sequences that are entirely composed of tandem repeats

## Pre-requisites

This workflow comes with a conda configuration file that will install most of its reprequisites. To use it, you must first install Anaconda or miniconda. I recommend miniconda. Download the appropriate installer from [here](https://docs.conda.io/en/latest/miniconda.html). If you are working on an MacOS device with an M1 or M2 chip, see the note below.

You will also need git, which is standard on most systems. Alternatively, it can be installed using conda:

    conda install git

## Installation

Clone this repository to your local system:

    git clone https://github.com/jmeppley/concatemer_finding

Navigate to the repository's root dir:

    cd concatemer_finding

Create a conda envirnoment using the included configuration:

    conda env create -p ./conda.env -f conda.yaml

Activate the environment

    conda activate ./conda.env

Just run snakemake to test. (-j 5: use 5 threads)

    snakemake -j 5
    
### A Note on conda for Mac M1 and M2 chips

Apple MacOS devices with an M1 or M2 processor can run either the "macOS Apple M1" (aka "arm") version of conda or the "macOS Intel x86" (aka "intel") version. The binaries compiled for intel will run a little more slowly than the M1 binaries, but not everything is avaiable for M1/M2 yet. In general, the bioconda packages are not yet available for M1/M2, but for this project, the one we need (minimap) is available in my conda channel. 

TL;DR: Either of the mac-os 64-bit versions of miniconda should work here.

## Run

Change the default configuration for snakemake. This assumes you run the repo
dir. See the snakemake docs for running elsewhere.

    snakemake --config output_dir='/path/to/output_dir' \
                       reads_template='/path/to/{sample}.reads.fasta'

The above assumes reads from each sample are in a single fasta file per sample.
