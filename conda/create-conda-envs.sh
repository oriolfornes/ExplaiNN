#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Snakemake
conda create -n snakemake -c bioconda snakemake
conda activate snakemake

# RSAT-core
conda create -n rsat -c bioconda -c conda-forge rsat-core=2020.02.29 \
    python=3.6.11

# Perl TFBS module
conda create -n tfbs -c bioconda -c conda-forge -c eumetsat \
    perl=5.26.2=h470a237_0 perl-bioperl-core perl-pdl
conda activate tfbs
git clone https://github.com/ComputationalRegulatoryGenomicsICL/TFBS.git
cd TFBS/
perl Makefile.PL && make && make install
cd ..
conda deactivate
