#!/usr/bin/env bash

SRC_DIR=${HOME}/motif-clustering
CHIP_SEQ_DIR=${HOME}/CAM/results/IRF4/CAM/ChIP-seq

for LABEL in "T95R" "WT"; do
    ${SRC_DIR}/meme2clusters.py \
        --out-dir ./ChIP-seq.${LABEL} \
        ${CHIP_SEQ_DIR}/${LABEL}/motifs/filters.meme
done
