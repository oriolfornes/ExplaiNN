#!/usr/bin/env bash

grep 2_PAD ~/CANN/resources/IRF4/ChIP-seq/01CLImagine_IRF4_intervals_noPAD156.csv | \
	awk 'BEGIN { FS=";" } {print "chr"$3"\t"$4"\t"$5}' | bedtools sort > WT.bed
grep 4_PAD ~/CANN/resources/IRF4/ChIP-seq/01CLImagine_IRF4_intervals_noPAD156.csv | \
	awk 'BEGIN { FS=";" } {print "chr"$3"\t"$4"\t"$5}' | bedtools sort > T95R.bed
grep 2_PAD ~/CANN/resources/IRF4/ChIP-seq/01CLImagine_IRF4_intervals_noPAD156.csv | \
        awk 'BEGIN { FS=";" } {print "chr"$3"\t"$4+$6-1-100"\t"$4+$6+100}' | bedtools sort > WT.201bp.bed
grep 4_PAD ~/CANN/resources/IRF4/ChIP-seq/01CLImagine_IRF4_intervals_noPAD156.csv | \
        awk 'BEGIN { FS=";" } {print "chr"$3"\t"$4+$6-1-100"\t"$4+$6+100}' | bedtools sort > T95R.201bp.bed
