#!/bin/bash

#SBATCH --nodes=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/DMRG-%A.out
#SBATCH --error=logs/DMRG-%A.err



exec=${1}
inputfile=${2}
inputbase=$(basename $inputfile .cfg)
nmin=${3}
nmax=${4}
outdir=logs/seed_$nmin-$nmax
mkdir -p $outdir

echo "Running job $SLURM_JOB_ID, seeds [$nmin - $nmax] at $HOSTNAME with inputfile $inputfile"
echo "CPUS ON  NODE  : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE  : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK  : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU    : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE   : $SLURM_MEM_PER_NODE"


parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log $exec $inputfile 1811 {} ::: $(seq $nmin $nmax)
#parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log $exec $inputfile {} ">" $outdir/$inputbase_{}.out ::: $(seq $nmin $nmax)


# parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log  $exec $inputfile {} '2>&1' ">" $outdir/$inputbase_{}.out ::: $(seq $nmin $nmax)

# '2>&1' ">" $outdir/$inputbase_{}.out
# Explanation

# $(seq $nmin $nmax) will print a sequence of seeds from nmin to nmax. Each number is piped with | to parallel,
# where the number is found in {1}.
# The lines are piped with "|" into GNU parallel, which takes a list as input argument. Each element of the list
# is expanded with {}.

# For each line, parallel will run exec with that line as input argument. The output of each execution is piped with ">"
# into file $outdir/$inputbase_{}.out   which will expand to something like  logs/jobarray_0/mbl_blabla_123.out

# There are many, many more arguments you can give to parallel, to control the number of cores, memory requirements per thread etc.