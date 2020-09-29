#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environments
conda create -n JASPAR-MoDisco -c bioconda -c conda-forge biopython=1.78  \
    dash=1.16.2 logomaker=0.8 plotly=4.10.0 python=3.6.11 rsat-core=2020.02.29 \
    snakemake=4.0.0
conda create -n meme -c bioconda/label/cf201901 -c conda-forge meme=5.0.2 \
    icu=58.2 libiconv=1.16

# Copy centrimo binaries to scripts folder
conda activate meme
BIN="$( dirname $( which centrimo ) )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# for FILE in "centrimo" "jaspar2meme"
for FILE in "centrimo"
do
    if [ -L $DIR/../workflow/scripts/centrimo ]; then
        rm $DIR/../workflow/scripts/centrimo
    fi
    ln -s $BIN/centrimo $DIR/../workflow/scripts/centrimo
done

# # Install PWMScan binaries
# # https://ccg.epfl.ch/pwmtools/pwmscan.php
# conda activate JASPAR-MoDisco
# # From https://stackoverflow.com/questions/59895/
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# cd $DIR/../workflow/scripts/pwmscan
# mkdir -p bin
# make clean && make cleanbin
# make && make install