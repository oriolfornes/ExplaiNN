#!/usr/bin/env bash

SRC_DIR=${HOME}/motif-clustering
CHIP_SEQ_DIR=${HOME}/CAM/results/IRF4/CAM/ChIP-seq

for LABEL in "T95R" "WT"; do
    for EXPERIMENT in 1; do
        ${SRC_DIR}/meme2clusters.py \
            --out-dir ./ChIP-seq/${LABEL}.${EXPERIMENT} \
            ${CHIP_SEQ_DIR}/${LABEL}.${EXPERIMENT}/motifs/filters.meme
    done
done

LABEL="T95R.1+WT.1"
${SRC_DIR}/meme2clusters.py --out-dir ./ChIP-seq/${LABEL} ${LABEL}.meme
