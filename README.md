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

Apple MacOS devices with an M1 or M2 processor can run either the "macOS Apple M1" (aka "arm") version of conda or the "macOS Intel x86" (aka "intel") version. The binaries compiled for intel will run more slowly than the M1 binaries, but not everything is avaiable for M1/M2 yet. For this project, the bioconda tools minimap2 and snakemake are missing. 

To use this repository on Apple Silicon (aka M1/M2 chips), you have two choices:

 1) Install the Intel conda from the start. This is the simplest option, but things may run a little slower.

 2) Install the Apple M1 conda, but configure it to use the Intel binaries when you need something specific. AFter you have installed conda, when you need to install something that's not available for M1/M2, you can configure conda to also look for intel binaries with:subdirectories for binaries:

````
conda config --add subdirs osx-64
````

I would recommend reverting this change when you're done here:

    conda config --remove subdirs osx-64

It's best to use primarily the arm binaries and only fall back to the intel ones when absolutely necessary. Avoid mixing binary types in your promary conda environment(s). Only try the intel binaries in separate, dedicated-use environments. 

## Run

Change the default configuration for snakemake. This assumes you run the repo
dir. See the snakemake docs for running elsewhere.

    snakemake --config output_dir='/path/to/output_dir' \
                       reads_template='/path/to/{sample}.reads.fasta'

The above assumes reads from each sample are in a single fasta file per sample.
