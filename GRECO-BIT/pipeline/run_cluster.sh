# variables
SNAKEFILE='Snakefile'
CONFIGFILE='config.yaml'

LOGDIR='/home/rvdlee/jaspar_2020/motif_discovery_pipeline/ModERN/logs/'
JOBNAME=RSAT

mkdir -p $LOGDIR
source activate py3

# PART 1
snakemake --cluster 'qsub -V -w e -N $JOBNAME -l h_vmem=2G' --snakefile $SNAKEFILE --configfile $CONFIGFILE --forcerun --drmaa-log-dir $LOGDIR
mv ~/$JOBNAME.o* $LOGDIR
mv ~/$JOBNAME.e* $LOGDIR

# # PART 2
# snakemake --cluster 'qsub -V -w e -N $JOBNAME -l h_vmem=2G' --snakefile $SNAKEFILE --configfile $CONFIGFILE
# mv ~/$JOBNAME.o* $LOGDIR
# mv ~/$JOBNAME.e* $LOGDIR

## RUN # qsub *sh
## MONITOR
# qstat -f
# watch -n 0.2 qstat -f
# logs are in home folder
# head RSAT.o* | less
# head RSAT.e* | less
# snakemake --cluster-status 10122 # <JOB_ID>
