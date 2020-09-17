#!/bin/bash

if [ "$#" -ne 1 ];
then
    echo
    echo "Usage: ./get_genome.sh genome"
    echo
else
    URL="https://hgdownload.soe.ucsc.edu/goldenPath/${1}/bigZips"
    echo
    echo "*** Downloading ${1}..."
    echo
    curl -O ${URL}/$1.fa.gz
    gunzip $1.fa.gz
    curl -O ${URL}/$1.chrom.sizes
fi
