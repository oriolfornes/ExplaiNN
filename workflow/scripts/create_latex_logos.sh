#!/usr/bin/env bash

################################################################################
#
################################################################################

usage="""
\n
$0 -l <latex header file> -i <input file> -o <output file>
\n
""";

\usepackage[export]{adjustbox}

TEMP=`getopt -o ho:i:l: -- "$@"`;

if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi
eval set -- "$TEMP";

input="NA";
output="NA";
latexheader="NA";
while true
do
  case "$1" in
    -i) input=$2; shift 2;;
    -o) output=$2; shift 2;;
    -l) latexheader=$2; shift 2;;
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

if [ ! -e $latexheader ]
then
    echo -e "\nPlease provide the latex header file\n";
    exit 1;
fi;

#BIN=/storage/mathelierarea/processed/jamondra/Projects/JASPAR/post_processing/jaspar_2020/motif_discovery_pipeline/bin
#cat $BIN/latex_header.txt > $output;
#latexheader="/home/rvdlee/JASPAR/jaimicore-jaspar_2020-a8db6feb9e0e/motif_discovery_pipeline/bin/latex_header.txt";
cat $latexheader > $output;

more $input | \
while read line
do
    tfname=$(echo $line | cut -d ' ' -f 2 | tr '_' '-');
    exp_ID=$(echo $line | cut -d ' ' -f 9 | tr '_' '-');
    #centrimo_file=$(echo $line | cut -d ' ' -f 16);
    centrimo_pval=$(echo $line | cut -d ' ' -f 10);
    
    motif_logo_original=$(echo $line | cut -d ' ' -f 8);
    # echo "===========>  $motif_logo_original";
    motif_logo_original_dirname=$(dirname $motif_logo_original);
    # echo "===========>  $motif_logo_original_dirname";
    # prevent error: ! LaTeX Error: Unknown graphics extension: .16_LE_WA_peak-motifs_m5_logo.png.
    # This removes any internal dots in the basename of the file, except for the dot in the extension (e.g. except the .png part)
    motif_logo_new=$(basename `echo $motif_logo_original` | perl -pe 's/(\.)(?=[^.]*\.)/\_/');
    # echo "===========>  $motif_logo_new";
    motif_logo=${motif_logo_original_dirname}/${motif_logo_new}
    # echo "===========>  $motif_logo";
    # this create a symbolic link so that that new motif filename will actually exist
    ln -s $motif_logo_original $motif_logo

    motif_pdf=$(echo $line | cut -d ' ' -f 11);
    
    # tfname=$(echo $line | cut -f3 | tr '_' '-');
    # exp_ID=$(echo $line | cut -f1);
    # centrimo_file=$(echo $line | cut -f4);
    #centrimo_pval=$(echo $line | cut -f5);
    # motif_logo=$(echo $line | cut -f6);
    # motif_pdf=$(echo $line | cut -f7);

   # if (( $(echo "$centrimo_pval < 0." | bc -l) ))
    #then
        #echo "TFname "$tfname;
        #echo "file "$centrimo_file;
        #echo "pval "$centrimo_pval;
    echo "\\section*{$tfname}";
    echo "\\section*{$exp_ID}";
    echo "Centrimo p-value = $centrimo_pval \\\\";
    #echo "Motif logo = $motif_logo \\\\";
        # pwmdir=$(dirname $centrimo_file);
        # echo -n "\\path{$pwmdir} \\\\";
        # pwmname=$(basename $centrimo_file .pssm.501bp.fa.sites.centrimo);
    echo "\\includegraphics{$motif_logo}";
    #fi;
done >> $output;
echo "\\end{document}" >> $output;
outdir=$(dirname $output);
pdflatex -output-directory $outdir $output;
