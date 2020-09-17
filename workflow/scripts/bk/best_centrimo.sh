#!/usr/bin/env bash

################################################################################
#
################################################################################

usage="""
\n
$0 -i <rsat dir>
\n
""";

TEMP=`getopt -o hi: -- "$@"`;
if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi
eval set -- "$TEMP";

rsat_dir=""
while true
do
  case "$1" in
    -i) rsat_dir=$2; shift 2;;
    -h) echo -e $usage; exit;;
    --) shift; break;;
    *) echo "Internal error!"; exit 1;;
  esac
done

# Checking options and parameters

#for centri_file in $rsat_dir/results/discovered_motifs/*_*nt*/peak*.centrimo
for centri_file in $rsat_dir/*.centrimo
do
    echo -e -n "$centri_file\t";
    cat $centri_file;
done > centri_$PPID;
awk 'BEGIN{mini=1} {if($2 < mini){mini=$2; minifile=$1;}} END{print minifile"\t"mini}' \
    centri_$PPID;
rm centri_$PPID;
