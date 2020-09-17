#!/usr/bin/env bash

################################################################################
#
################################################################################

usage="""
\n
$0 -r <repository>
\n
""";

TEMP=`getopt -o r:h -- "$@"`;
if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi;
eval set -- "$TEMP";

repo="";
while true
do
  case "$1" in
    -r) repo=$2; shift 2;; 
    -h) echo -e $usage; exit;;
    --) shift; break;;
    *) echo "Internal error!"; exit 1;;
  esac;
done;

# Checking options and parameters

if [ ! -d $repo ]
then
    echo "##### Repository $repo does not exist! #####";
    exit;
fi;

# MAIN

for peakfile in $repo/*1001bp.fa
do
    tf=$(basename $peakfile .1001bp.fa);
    centrimo_pval_motif1=$(cat $peakfile.homer_motif1.best_hits.tsv.centrimo_pval.txt);
    centrimo_pval_motif1=$(printf "%f" $centrimo_pval_motif1);
    centrimo_pval_motif2=$(cat $peakfile.homer_motif2.best_hits.tsv.centrimo_pval.txt);
    centrimo_pval_motif2=$(printf "%f" $centrimo_pval_motif2);
    centrimo_pval_motif3=$(cat $peakfile.homer_motif3.best_hits.tsv.centrimo_pval.txt);
    centrimo_pval_motif3=$(printf "%f" $centrimo_pval_motif3);
    best_motif="NA";
    best_centrimo="NA";
    if (( $(echo "$centrimo_pval_motif1 < $centrimo_pval_motif2" | bc -l) ))
    then
        if (( $(echo "$centrimo_pval_motif1 < $centrimo_pval_motif3" | bc -l) ))
        then
            if (( $(echo "$centrimo_pval_motif1 < 0" | bc -l) ))
            then
                best_motif=$repo/$tf.101bp.fa.homer/homerResults/motif1.hits;
                best_centrimo=$centrimo_pval_motif1;
            fi;
        else
            if (( $(echo "$centrimo_pval_motif3 < 0" | bc -l) ))
            then
                best_motif=$repo/$tf.101bp.fa.homer/homerResults/motif3.hits;
                best_centrimo=$centrimo_pval_motif3;
            fi;
        fi;
    else
        if (( $(echo "$centrimo_pval_motif2 < $centrimo_pval_motif3" | bc -l) ))
        then
            if (( $(echo "$centrimo_pval_motif2 < 0" | bc -l) ))
            then
                best_motif=$repo/$tf.101bp.fa.homer/homerResults/motif2.hits;
                best_centrimo=$centrimo_pval_motif2;
            fi;
        else
            if (( $(echo "$centrimo_pval_motif3 < 0" | bc -l) ))
            then
                best_motif=$repo/$tf.101bp.fa.homer/homerResults/motif3.hits;
                best_centrimo=$centrimo_pval_motif3;
            fi;
        fi;
    fi;
    if [ $best_motif != "NA" ]
    then
        /lsc/weblogo/seqlogo -f $best_motif.fa -o $best_motif -x position -y bits \
            -c -p -n -Y -F PDF;
        echo -e "$(basename $tf _peaks.narrowPeak)\t$best_centrimo\t$best_motif.pcm\t$best_motif.pdf";
    else
        echo -e "$(basename $tf _peaks.narrowPeak)\tNA\tNA\tNA";
    fi;
done;
