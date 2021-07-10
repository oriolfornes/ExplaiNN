#!/usr/bin/env bash

# Initialize
CNN_UNITS=16
MOTIF_LENGTH=24
ROOT_DIR=${HOME}/Work/deep-motif-discovery
CAM_DIR=${ROOT_DIR}/workflow/CAM
OUTPUT_DIR="${ROOT_DIR}/results/DREAM5/CAM.cnn-units=${CNN_UNITS}&motif-length=${MOTIF_LENGTH}"
SEQUENCES_DIR=${ROOT_DIR}/results/DREAM5/FASTA

# CAM
eval "$(conda shell.bash hook)"
source activate deep-motif-discovery-pipeline

cd ${CAM_DIR}

for F in `ls ${SEQUENCES_DIR}/Train`; do

    TF=${F: 0: -6}

    if ! [[ -d ${OUTPUT_DIR}/${TF} ]]; then

        ./train.py \
            --output-dir ${OUTPUT_DIR}/${TF} \
            --cnn-units ${CNN_UNITS} \
            --input-data linear \
            --motif-length ${MOTIF_LENGTH} \
            ${SEQUENCES_DIR}/Train/${TF}.fa.gz \
            ${SEQUENCES_DIR}/Test/${TF}.fa.gz

    fi

done
