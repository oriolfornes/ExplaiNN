#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n JASPAR-MoDisco -c bioconda -c conda-forge biopython logomaker \
    python=3.6.11 rsat-core=2020.02.29 snakemake=4.0.0

# Install PWMScan


