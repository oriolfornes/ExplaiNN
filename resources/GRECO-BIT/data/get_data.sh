#!/bin/bash

INTERVALS_DIR="/home_local/vorontsovie/greco-data/release_3.2020-08-08/chipseq/results/train_intervals/"
INTERVAL_FILES=$(ls ${INTERVALS_DIR})
FILE_EXTENSION=".chipseq.train.interval"

for INTERVAL_FILE in ${INTERVAL_FILES}; do
    INTERVAL_DIR=${INTERVAL_FILE:0:${#INTERVAL_FILE}-${#FILE_EXTENSION}}
    mkdir ${INTERVAL_DIR}
    ln -s ${INTERVALS_DIR}/${INTERVAL_FILE} ${INTERVAL_DIR}
done