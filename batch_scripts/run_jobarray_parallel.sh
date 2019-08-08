#!/bin/bash

#SBATCH --nodes=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/DMRG-%A_%a.out
#SBATCH --error=logs/DMRG-%A_%a.err



exec=${1}
inputfile=${2}
inputbase=$(basename $inputfile .cfg)
nmin=$SLURM_ARRAY_TASK_MIN
nmax=$((nmin + SLURM_ARRAY_TASK_STEP - 1))
outdir=logs/jobarray_$SLURM_ARRAY_JOB_ID
mkdir -p $outdir

echo "Running job $SLURM_ARRAY_JOB_ID, id [$nmin - $nmax] at $HOSTNAME with inputfile $inputfile"
echo "ARRAY JOB ID     : $SLURM_ARRAY_JOB_ID"
echo "ARRAY TASK COUNT : $SLURM_ARRAY_TASK_COUNT"
echo "ARRAY TASK ID    : $SLURM_ARRAY_TASK_ID"
echo "ARRAY TASK MAX   : $SLURM_ARRAY_TASK_MAX"
echo "ARRAY TASK MIN   : $SLURM_ARRAY_TASK_MIN"
echo "ARRAY TASK STEP  : $SLURM_ARRAY_TASK_STEP"

echo "CPUS ON  NODE  : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE  : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK  : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU    : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE   : $SLURM_MEM_PER_NODE"


parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log --results $outdir $exec $inputfile {} ::: $(seq $nmin $nmax)

#parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log  $exec $inputfile {} '2>&1' ">" $outdir/$inputbase_{}.out ::: $(seq $nmin $nmax)

# '2>&1' ">" $outdir/$inputbase_{}.out
# Explanation

# $(seq $nmin $nmax) will print a sequence of seeds from nmin to nmax. Each number is piped with | to parallel,
# where the number is found in {1}.
# The lines are piped with "|" into GNU parallel, which takes a list as input argument. Each element of the list
# is expanded with {}.

# For each line, parallel will run exec with that line as input argument. The output of each execution is piped with ">"
# into file $outdir/$inputbase_{}.out   which will expand to something like  logs/jobarray_0/mbl_blabla_123.out

# There are many, many more arguments you can give to parallel, to control the number of cores, memory requirements per thread etc.