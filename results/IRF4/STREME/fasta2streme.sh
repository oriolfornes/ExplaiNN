#!/usr/bin/env bash

#streme --p WT.fa --oc WT --maxw 21
#streme --p T95R.fa --oc T95R --maxw 21 
#streme --p WT.201bp.fa --oc WT.201bp --maxw 21
#streme --p T95R.201bp.fa --oc T95R.201bp --maxw 21
streme --p WT-unique.fa --oc WT-unique --maxw 21
streme --p T95R-unique.fa --oc T95R-unique --maxw 21
streme --p T95R-WT-intersect.fa --oc T95R-WT-intersect --maxw 21
