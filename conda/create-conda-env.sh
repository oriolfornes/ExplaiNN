#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n JASPAR-MoDisco -c bioconda -c conda-forge biopython=1.78 \
    logomaker=0.8 python=3.6.11 rsat-core=2020.02.29 snakemake=4.0.0

# Install PWMScan binaries
# https://ccg.epfl.ch/pwmtools/pwmscan.php
conda activate JASPAR-MoDisco
# From https://stackoverflow.com/questions/59895/
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR/../workflow/scripts/pwmscan
make && make install