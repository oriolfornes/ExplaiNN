#!/bin/bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Extract peak summits: chromosome, start, end
# In the GRECO-BIT chipseq.train.interval files, column 1 corresponds to chrom
# and column 4 to the summit position.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DIR="/home_local/vorontsovie/greco-data/release_3.2020-08-08/chipseq/results/train_intervals"
FILES=$(ls ${DIR})
EXT=".chipseq.train.interval"

for FILE in ${FILES}; do
    NEW_DIR=${FILE:0:${#FILE}-${#EXT}}
    mkdir ${NEW_DIR}
    NEW_FILE="${NEW_DIR}/${NEW_DIR}_peak_summits.bed"
    cut -f 1,4 ${DIR}/${FILE} | \
    grep -v "^#" | \
    awk '{print($1"\t"$2-1"\t"$2);}' | \
    LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ${NEW_FILE}
done