#!/bin/bash

SNAKE_FILE="Snakefile"
CONFIG_FILE="config.yaml"

snakemake --cores 4 --snakefile ${SNAKE_FILE} --configfile ${CONFIG_FILE}
