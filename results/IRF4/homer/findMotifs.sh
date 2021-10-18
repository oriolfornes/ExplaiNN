#!/usr/bin/env bash

findMotifs.pl WT.201bp.fa fasta WT.201bp.homer-findMotifs -fasta WT.201bp.scrambled.fa
findMotifs.pl T95R.201bp.fa fasta T95R.201bp.homer-findMotifs -fasta T95R.201bp.scrambled.fa
findMotifs.pl WT-unique.fa fasta WT-unique.homer-findMotifs -fasta WT-unique.scrambled.fa
findMotifs.pl T95R-unique.fa fasta T95R-unique.homer-findMotifs -fasta T95R-unique.scrambled.fa
findMotifs.pl T95R-WT-intersect.fa fasta T95R-WT-intersect.homer-findMotifs -fasta T95R-WT-intersect.scrambled.fa
