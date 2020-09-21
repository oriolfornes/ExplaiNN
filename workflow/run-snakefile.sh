#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate JASPAR-MoDisco

# Run workflow
snakemake --snakefile Snakefile --cores 16 --configfile ../config/config.yml