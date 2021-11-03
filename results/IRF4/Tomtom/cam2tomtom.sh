#!/usr/bin/env bash

CAM_DIR=${HOME}/CAM/results/IRF4/CAM

for ASSAY in "Affi-seq" "ChIP-seq" "HT-SELEX"; do
    for LABEL in "T95R" "WT"; do
        tomtom -no-ssc \
            -verbosity 1 \
            -min-overlap 5 \
            -dist pearson \
            -evalue \
            -thresh 10.0 \
            -oc CAM.${ASSAY}.${LABEL} \
            ${CAM_DIR}/${ASSAY}/${LABEL}/motifs/filters.meme \
            JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
    done
done

CAM_DIR=${HOME}/CAM/results/IRF4/CAM/ChIP-seq
WT=${CAM_DIR}/WT/motifs/filters.meme
T95R=${CAM_DIR}/T95R/motifs/filters.meme
tomtom -no-ssc \
    -verbosity 1 \
    -min-overlap 5 \
    -dist pearson \
    -evalue \
    -thresh 10.0 \
    -oc CAM.ChIP-seq.T95R-WT \
    ${WT} ${T95R}
