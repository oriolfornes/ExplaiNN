#!/usr/bin/env bash

BASE_DIR=${HOME}/CAM

fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA1_A07_AC40NGTATGA_Read1.fastq.gz -o ./HT-SELEX/YWIA1_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA2_A07_AC40NGTATGA_Read1.fastq.gz -o ./HT-SELEX/YWIA2_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA3_A07_AC40NGTATGA_Read1.fastq.gz -o ./HT-SELEX/YWIA3_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA1_B07_AT40NGAGAGG_Read1.fastq.gz -o ./HT-SELEX/YWIA1_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA2_B07_AT40NGAGAGG_Read1.fastq.gz -o ./HT-SELEX/YWIA2_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/HT-SELEX/YWIA3_B07_AT40NGAGAGG_Read1.fastq.gz -o ./HT-SELEX/YWIA3_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB1_A07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB1_A07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB1_A07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB1_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB2_A07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB2_A07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB2_A07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB2_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB3_A07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB3_A07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB3_A07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB3_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB1_B07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB1_B07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB1_B07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB1_B07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB2_B07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB2_B07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB2_B07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB2_B07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB3_B07_AffiSeq_Read1.fastq.gz -o ./Affi-seq/YWIB3_B07_AffiSeq_Read1.fastq.gz -I ${BASE_DIR}/resources/IRF4/Affi-seq/YWIB3_B07_AffiSeq_Read2.fastq.gz -O ./Affi-seq/YWIB3_B07_AffiSeq_Read2.fastq.gz -A -G -w 8
