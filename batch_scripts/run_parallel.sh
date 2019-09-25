#!/bin/bash

#SBATCH --nodes=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err




PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-e <executable>                     : Path to executable (default = "")
-f <simfile>                        : Path to simulation file, two columns formatted as [configfile seed] (default = "")
-o <output logfile>                 : Path to GNU parallel logfile (default = "")
EOF
  exit 1
}

while getopts he:f:o: o; do
    case $o in
        (h) usage ;;
        (e) exec=$OPTARG;;
        (f) simfile=$OPTARG;;
        (o) outfile=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


if [ ! -f $simfile ]; then
    echo "Simfile is not a valid file: $simfile"
    exit 1
else

if [ -z "$outfile" ] ; then
    simbase=$(basename $simfile .sim)
    outfile=logs/$simbase.log
fi



echo "Running job $SLURM_JOB_ID at $HOSTNAME with simfile $simfile"
outdir=$(dirname $outfile)
mkdir -p $outdir


num_cols=$(awk '{print NF}' $simfile | head -n 1)

if [ "$num_cols" -eq 2 ]; then
    cat $simfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outfile --colsep ' ' "$exec -i {1} -r {2} &> $outdir/${simbase}/{1}_{2}.out"
elif [ "$num_cols" -eq 3 ]; then
    cat $simfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outfile --colsep ' ' "$exec -i {1} -r {2} -s {3} &> $outdir/${simbase}/{1}_{2}_{3}.out"
else
    echo "Case not implemented"
    exit 1
fi


#cat $simfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outfile "$exec -i $inputfile -r {} &> $outdir/${inputbase}_{}.out"

#
#inputbase=$(basename $inputfile .cfg)
#mkdir -p $outdir
#
#echo "CPUS ON  NODE  : $SLURM_CPUS_ON_NODE"
#echo "CPUS PER NODE  : $SLURM_JOB_CPUS_PER_NODE"
#echo "CPUS PER TASK  : $SLURM_CPUS_PER_TASK"
#echo "MEM PER CPU    : $SLURM_MEM_PER_CPU"
#echo "MEM PER NODE   : $SLURM_MEM_PER_NODE"
#
#if [ -n "$seedfile" ]; then
#    seedbase=$(basename $seedfile .txt)
#    echo "Running job $SLURM_JOB_ID, with seeds from file $seedfile at $HOSTNAME with inputfile $inputfile"
#    cat $seedfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase_$seedbase.log "$exec -i $inputfile -r {} &> $outdir/${inputbase}_{}.out"
#elif [ -n "$nmin" ] && [ -n "$nmax" ] ;then
#    echo "Running job $SLURM_JOB_ID, seeds [$nmin - $nmax] at $HOSTNAME with inputfile $inputfile"
#    parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log "$exec -i $inputfile -r {} &> $outdir/${inputbase}_{}.out" ::: $(seq $nmin $nmax)
#fi
#

# For initial state enumeration
#realization=1811
#parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log "$exec -i $inputfile -r $realization -s {} &> $outdir/${inputbase}_${realization}_{}.out" ::: $(seq $nmin $nmax)






# parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log  $exec $inputfile {} '2>&1' ">" $outdir/$inputbase_{}.out ::: $(seq $nmin $nmax)

# Explanation

# $(seq $nmin $nmax) will print a sequence of seeds from nmin to nmax. Each number is piped with | to parallel,
# where the number is found in {1}.
# The lines are piped with "|" into GNU parallel, which takes a list as input argument. Each element of the list
# is expanded with {}.

# For each line, parallel will run exec with that line as input argument. The output of each execution is piped with ">"
# into file $outdir/$inputbase_{}.out   which will expand to something like  logs/jobarray_0/mbl_blabla_123.out

# There are many, many more arguments you can give to parallel, to control the number of cores, memory requirements per thread etc.