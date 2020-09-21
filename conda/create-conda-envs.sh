#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n JASPAR-MoDisco -c bioconda -c conda-forge -c Eumetsat \
    perl=5.26.2=h470a237_0 perl-pdl=2.019 python=3.6.11 rsat-core=2020.02.29 \
    snakemake=4.0.0

# Install TFBS module
conda activate JASPAR-MoDisco
git clone https://github.com/ComputationalRegulatoryGenomicsICL/TFBS.git
cd TFBS/
perl Makefile.PL && make && make install
cd ..
conda deactivate
