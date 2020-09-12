# JASPAR 2020

Depending on the ChIP-seq data repository (e.g., ChIP-atlas or ENCODE) the pipeline may vary slightly. For this reason there is a folder for each repository with its respective Snakefile and Config file.

This folder is divided in 2 subfolders:

- bin       : with all the executable files called in the workflow.
- ChIP-atlas: contains all the analysis performed in ChIP-atlas datasets + Snakefile + Config file
- ModERN: contains the files for processing the ModERN worm and fly datasets described in http://www.genetics.org/content/208/3/937


## Contributors
Author: Jaime A Castro-Mondragon  
Contributors: Robin van der Lee  


## Dependencies.

The following software and R libraries must be installed in order to run the complete pipeline.

- *R*
	- data.table
	- dplyr
	- future.apply
	- reshape2  
```install.packages("data.table"); install.packages("dplyr"); install.packages("future.apply"); install.packages("reshape2")```  

- *RSAT* (Regulatory Sequences analysis Tools) software: http://rsat-tagc.univ-mrs.fr/rsat/
- Perl TFBS modules: http://tfbs.genereg.net/
- *snakemake*: https://snakemake.readthedocs.io/en/stable/index.html
- pdfunite

## Download ChIP-seq datasets from ChIP-atlas

1. Change the directory:
```
cd motif_discovery_pipeline/ChIP-atlas
```


2. Launch the script. You must specify the directory (current directory '.' in this example) and the genome version (ce10, C elegans). This program will first download the experiment table from ChIP-atlas and then all the experiments (narrowpeak format) for the given genome. This script checks the file size after the download, empty files are removed.
```
Rscript bin/ChIP_atlas_download_TF_peaks_by_genome.R . ce10     ## C elegans
Rscript bin/ChIP_atlas_download_TF_peaks_by_genome.R . dm3      ## D melanogaster
Rscript bin/ChIP_atlas_download_TF_peaks_by_genome.R . sacCer3  ## S Cereviseae
```


## Download genomic files: chromosome sizes and fasta files. C elegans example.


1. Go to UCSC download page: http://hgdownload.cse.ucsc.edu/downloads.html#c_elegans

2. Look for the ce10 (or any genome) table 'Oct. 2010 (WS220/ce10)', click on 'Full data set'

3. Download the genome size file:
   ```
   wget http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.chrom.sizes
   ```

4. Download the ce10 genome in 2bit format:
   ```
   wget http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit
   ```

5. Download the script 2bittofa, this script convert the 2bit file into fasta sequences:
   ```
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
   ```

6. Convert the 2bit file to fasta sequences:
   ```
   twoBitToFa ce10.2bit ce10.fa
   ```

7. twoBitToFa information: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/


```
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
chmod +x twoBitToFa
./twoBitToFa sacCer3.2bit sacCer3.fa
```

#### Genomes for the ModERN data 
```
http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.2bit
./twoBitToFa dm6.2bit dm6.fa
```

```
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS245.genomic.fa.gz
gunzip c_elegans.PRJNA13758.WS245.genomic.fa.gz 

# https://www.biostars.org/p/173963/ =>  pip install pyfaidx
faidx c_elegans.PRJNA13758.WS245.genomic.fa -i chromsizes > c_elegans.PRJNA13758.WS245.genomic.chrom.sizes

bash convert_chrlabels_c_elegans.sh 
```


## Configuration file

Set the following variables in the config file:


### Software variables

Assuming that you are in the *motif_discovery_pipeline/ChIP-atlas* folder:

- *bin* : ../bin 
- *python* : python2.7 
- *RSAT* : /lsc/rsat

You can see where RSAT program are using the following command: ```echo $RSAT```


### Organism specific variables

- data_folder : This is the folder where the input file provided as *narrowPeak* must be placed. See *Expected folder structure* for more information about the folder and file names.
- genome_fasta: The path to the *.fa* file with the organism's genome sequence. See *Download genomic files* section.
- genome_size: the path to the *.chrom.sizes* file. See *Download genomic files* section.
- TF_Experiment_map: The path to the experiment table. See *Download ChIP-seq datasets from ChIP-atlas*.
- out_dir: The path to results folder where all the analysis will be stored, one folder per experiment. Don't forget to indicate the genome in the path. Example: ``` ChIP-atlas_results/ce10```


## Launching the pipeline

### Expected folder structure

The Snakefile will launch a motif analysis + centrality enrichment test for each discovered motif of each dataset (experiment). To launch the pipeline we only required a *narrowPeak* file that must be located in the 'data_folder'specified in the config file (see *configuration file* section).

Example:
```
cd data_folder

.
├── ERX032305_S2_Lint-1
│   └── ERX032305_S2_Lint-1_peaks.narrowPeak
└── ERX032306_Kc_Lint-1
    └── ERX032306_Kc_Lint-1_peaks.narrowPeak
```

Every EXPERIMENT_FOLDER in the data folder must contain a narrowPeak file with the same name and the suffix *_peaks*, see tree above.
Example:
```
EXPERIMENT_FOLDER = ERX032305_S2_Lint-1
narrowPeak        = ERX032305_S2_Lint-1_peaks.narrowPeak
```
This is the data structure required by the pipeline.


### Launching the *snakemake* workflow

The pipeline must be ran in two steps:

- Step 1: motif discovery. The program *RSAT peak-motifs* is called, this program runs four motif discovery algorithms (over-represented k-mers, over-represented spaced-k-mers, positionally biased k-mers, k-mers enriched in central windows) to discover significant motifs, see 10.1093/nar/gkr1104. This can be launched by running the *Snakefile*, you can parallelize the process using the argument *--cores*.

The *Snakefile* contain a series of rules, for the motif discovery steps, it searches for datasets in *narrowPeak* format. 

Once completed, this step will generate one folder per dataset with the *RSAT peak-motifs* results. 

```
snakemake [--cores X]
```

- Step 2: selecting centrally enriched motifs. Note that given the ChIP-seq dataset quality and the number of peaks may vary among the datasets, many datasets may not produce any significant motif. A priori, the number of discovered motifs is unknow, to avoid the use of *dynamic files*, the pipeline will consider as input all the discovered motifs in *transfac* format resulting from the step 1,**following certain folder structure** (resulting from the rules ran in the step 1, see NOTE 1 in the *Snakefile*). For this reason, the pipeline must be ran in two steps.


```
snakemake [--cores X]
```


## Troubleshooting


### rule RSAT_peakmotifs_per_exp

- *RSAT peak-motifs* none motif found: verify that the input *NarrowPeak* is not empty. In the case of ChIP-atlas, many of the downloaded files are empty, to avoid this error, remove these files before launching the pipeline.

- *RSAT local-word-analysis* python version: this program is written in python and called within RSAT peak-motifs. In case that the default python environment is 3.5, local-word-analysis will not run, therefore the solution is modify directly the line in the $RSAT/perl-scripts/peak-motifs script adding "python2.7" before the program is called.

- *RSAT position-analysis* founds highly repetitive/uninformative motifs: The default interval for position-analysis is 50nt, which is too large for the JASPAR pipeline (where the peaks have a size of 101), this may produce non-relevant artifact motifs. By changing the interval to 25 the program may found the expected TF motif. Note that this change was already set in the config file.


### rule annotate_best_centrimo_experiment

- Verify that the variables in the *awk* command correspond to the correct fields in your experiment table. The current *awk* variables in the rule are specific for the ChIP-atlas experiment table.


### rule best_centrimo_experiment_logo

- Input files not found: this may occur when the Experiment_TF ID contains certain characters as: *(*, *)*, *[*, *]*, *.*, */*, " "(one or more spaces). Before launching the pipeline, verify that none of the folders in your data folder contains such characters, the simplest option is remove them or substitute them by '*-*'. Note that this step was already implemented for the ChIP-atlas datasets.


## Memory requirements

If you launch the *Snakefile* using multiple cores (See *Launching the *snakemake* workflow* section), each core will use ~1GB of memory in your computer.
