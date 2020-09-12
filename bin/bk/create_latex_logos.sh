#!/usr/bin/env bash

################################################################################
#
################################################################################

usage="""
\n
$0 -i <input file> -o <output file>
\n
""";

TEMP=`getopt -o ho:i: -- "$@"`;
if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi
eval set -- "$TEMP";

input="NA";
output="NA";
while true
do
  case "$1" in
    -i) input=$2; shift 2;;
    -o) output=$2; shift 2;;
    -h) echo -e $usage; exit;;
    --) shift; break;;
    *) echo "Internal error!"; exit 1;;
  esac
done

# Checking options and parameters

if [ ! -e $input ]
then
    echo -e "\nPlease provide the input file\n";
    exit 1;
fi;

BIN=/storage/mathelierarea/processed/jamondra/Projects/JASPAR/post_processing/bin
cat $BIN/latex_header.txt > $output;
sort -k3,3 -k5,5n $input | \
while read line
do
    tfname=$(echo $line | cut -d ' ' -f 3 | tr '_' '-');
    centrimo_file=$(echo $line | cut -d ' ' -f 4);
    centrimo_pval=$(echo $line | cut -d ' ' -f 5);
    motif_logo=$(echo $line | cut -d ' ' -f 6);
    if (( $(echo "$centrimo_pval < 0." | bc -l) ))
    then
        #echo "TFname "$tfname;
        #echo "file "$centrimo_file;
        #echo "pval "$centrimo_pval;
        echo "\\section*{$tfname}";
        echo "Centrimo p-value = $centrimo_pval \\\\";
	#echo "Motif logo = $motif_logo \\\\";
        # pwmdir=$(dirname $centrimo_file);
        # echo -n "\\path{$pwmdir} \\\\";
        # pwmname=$(basename $centrimo_file .pssm.501bp.fa.sites.centrimo);
        echo "\\includegraphics{$motif_logo}";
    fi;
done >> $output;
echo "\\end{document}" >> $output;
outdir=$(dirname $output);
pdflatex -output-directory $outdir $output;
