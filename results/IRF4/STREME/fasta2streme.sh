#!/usr/bin/env bash

FASTA_DIR=${HOME}/CAM/resources/IRF4/ChIP-seq

streme --p ${FASTA_DIR}/WT.fa --oc WT --maxw 21
streme --p ${FASTA_DIR}/T95R.fa --oc T95R --maxw 21
streme --p ${FASTA_DIR}/WT-unique.fa --oc WT-unique --maxw 21
streme --p ${FASTA_DIR}/T95R-unique.fa --oc T95R-unique --maxw 21
streme --p ${FASTA_DIR}/T95R-WT-intersect.fa --oc T95R-WT-intersect --maxw 21
