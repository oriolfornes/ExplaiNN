#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

#conda activate RSAT
conda activate rsat
snakemake --snakefile ./RSAT/Snakefile --cores 8 --configfile \
    ../config/config.yml