#!/usr/bin/env bash

bedtools getfasta -fi ~/CANN/resources/genomes/hg38/hg38.fa -bed WT.bed -fo WT.fa
bedtools getfasta -fi ~/CANN/resources/genomes/hg38/hg38.fa -bed T95R.bed -fo T95R.fa
bedtools getfasta -fi ~/CANN/resources/genomes/hg38/hg38.fa -bed WT.201bp.bed -fo WT.201bp.fa
bedtools getfasta -fi ~/CANN/resources/genomes/hg38/hg38.fa -bed T95R.201bp.bed -fo T95R.201bp.fa
