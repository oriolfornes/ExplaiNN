#!/usr/bin/env bash

fastp -i ../../resources/IRF4/FASTQ/YWIA1_A07_AC40NGTATGA_Read1.fastq.gz -o ./FASTQ/YWIA1_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIA2_A07_AC40NGTATGA_Read1.fastq.gz -o ./FASTQ/YWIA2_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIA3_A07_AC40NGTATGA_Read1.fastq.gz -o ./FASTQ/YWIA3_A07_AC40NGTATGA_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIA1_B07_AT40NGAGAGG_Read1.fastq.gz -o ./FASTQ/YWIA1_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIA2_B07_AT40NGAGAGG_Read1.fastq.gz -o ./FASTQ/YWIA2_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIA3_B07_AT40NGAGAGG_Read1.fastq.gz -o ./FASTQ/YWIA3_B07_AT40NGAGAGG_Read1.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB1_A07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB1_A07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB1_A07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB1_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB2_A07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB2_A07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB2_A07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB2_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB3_A07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB3_A07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB3_A07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB3_A07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB1_B07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB1_B07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB1_B07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB1_B07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB2_B07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB2_B07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB2_B07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB2_B07_AffiSeq_Read2.fastq.gz -A -G -w 8
fastp -i ../../resources/IRF4/FASTQ/YWIB3_B07_AffiSeq_Read1.fastq.gz -o ./FASTQ/YWIB3_B07_AffiSeq_Read1.fastq.gz -I ./resources/IRF4/FASTQ/YWIB3_B07_AffiSeq_Read2.fastq.gz -O ./FASTQ/YWIB3_B07_AffiSeq_Read2.fastq.gz -A -G -w 8