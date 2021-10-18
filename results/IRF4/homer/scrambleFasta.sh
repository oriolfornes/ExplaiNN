#!/usr/bin/env bash

scrambleFasta.pl WT.201bp.fa > WT.201bp.scrambled.fa
scrambleFasta.pl T95R.201bp.fa > T95R.201bp.scrambled.fa
scrambleFasta.pl WT-unique.fa > WT-unique.scrambled.fa
scrambleFasta.pl T95R-unique.fa > T95R-unique.scrambled.fa
scrambleFasta.pl T95R-WT-intersect.fa > T95R-WT-intersect.scrambled.fa
