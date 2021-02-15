#!/usr/bin/env bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environments
conda create -n jmodisco -c bioconda -c conda-forge -c pytorch        \
    biopython=1.78 click=7.1.2 click-option-group=0.5.1 dash=1.19.0   \
    genomepy=0.9.3 ignite=0.4.3 jupyterlab=3.0.7 logomaker=0.8        \
    pybedtools=0.8.1 python=3.6.12 seaborn=0.11.1 scikit-learn=0.24.1 \
    snakemake=5.3.0 tensorboard=2.4.1 torchvision=0.8.1 tqdm=4.56.0
conda create -n meme -c bioconda/label/cf201901 -c conda-forge        \
    meme=5.0.2 icu=58.2 libiconv=1.16

# Copy centrimo binaries to scripts folder
conda activate meme
BIN="$( dirname $( which meme ) )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
for FILE in "centrimo" "fasta-shuffle-letters"
do
    if [ -L $DIR/../workflow/scripts/${FILE} ]; then
        rm $DIR/../workflow/scripts/${FILE}
    fi
    ln -s $BIN/${FILE} $DIR/../workflow/scripts/${FILE}
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