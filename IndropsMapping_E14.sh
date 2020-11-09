logdir="filter_logs_redo/" # path for writing log files
yaml="MIAE14_022820.yaml" #Specify your specific yaml file with indices
nlibs=8 # number of libraries to process
nlanes=4 # number of lanes to process

######################################################
# In general, no need to change anything below here. #
######################################################
N1=$nlanes # number of workers for "filter" step
N3=$nlibs # number of workers for "sort" step
N4=200 # number of workers for "quantify" step

# bug recently introduced that messes up multiple
# workers in the aggregate step
# N5=$nlibs # number of workers for "aggregate" step

module load gcc/6.2.0
module load bowtie/1.2.2
module load java/jdk-1.8u112
module load rsem/1.3.0

mkdir -p ${logdir}

#Specify virtual python environment and indrops python script, then activate virtual python environment
pyenv=/n/groups/neuroduo/brian/pyndrops
indrops=/n/groups/neuroduo/brian/indrops/indrops.py
myPython=/n/groups/neuroduo/brian/pyndrops/bin/python
source ${pyenv}/bin/activate ${pyenv}


## PART 1 (filter)
# Use this if just one set of fastq files
jid1=$( sbatch -p medium --array 1-${N1} -n 3 -N 1 --job-name F --mem 30000 -t 1-00:00:00 -o "$logdir"filter.out -e "$logdir"filter.out --wrap """ ${myPython} ${indrops} $yaml filter --total-workers ${N1} --worker-index \$((\$SLURM_ARRAY_TASK_ID-1)) """ )
jid1=${jid1##* }
echo $jid1

## PART 2 (identify abundant barcodes)
jid2=$( sbatch --dependency=afterany:$jid1 -p short -n 1 -N 1 --job-name I --mem 60000 -t 4:00:00 -o "$logdir"abundant.out -e "$logdir"abundant.out --wrap """ ${myPython} $indrops $yaml identify_abundant_barcodes """ )
jid2=${jid2##* }
echo $jid2

## PART 3 (sort)
# max number of workers = number of libraries
jid3=$( sbatch --dependency=afterany:$jid2 -p short --array 1-${N3} -n 3 -N 1 --job-name S --mem 60000 -t 12:00:00 -o "$logdir"sort_worker_%a.out -e "$logdir"sort_worker_%a.out --wrap """ ${myPython} $indrops $yaml sort --total-workers ${N3} --worker-index \$((\$SLURM_ARRAY_TASK_ID-1)) """ )
jid3=${jid3##* }
echo $jid3

## PART 4 (quantify)
jid4=$( sbatch --dependency=afterany:$jid3 -p short --array 1-${N4} -n 3 -N 1 --job-name Q --mem 100000 -t 10:00:00 -o "$logdir"quant_worker_%a.out -e "$logdir"quant_worker_%a.out --wrap """ ${myPython} $indrops $yaml quantify --min-reads 250 --min-counts 0 --total-workers ${N4} --worker-index \$((\$SLURM_ARRAY_TASK_ID-1)) --no-bam """)
jid4=${jid4##* }
echo $jid4

## PART 5 (aggregate)
jid5=$( sbatch --dependency=afterany:$jid4 -p short -n 4 -N 1 --job-name A --mem 25000 -t 05:00:00 -o "$logdir"agg.out -e "$logdir"agg.out --wrap """ ${myPython} $indrops $yaml aggregate --total-workers ${N4} --no-bam """ )
jid5=${jid5##* }
echo $jid5
