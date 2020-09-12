#!/bin/bash

SNAKE_FILE="Snakefile"
CONFIG_FILE="config.yaml"

snakemake --cores 8 --snakefile ${SNAKE_FILE} --configfile ${CONFIG_FILE}
