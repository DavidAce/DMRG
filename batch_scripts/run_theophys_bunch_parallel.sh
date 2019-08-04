#!/bin/bash
#SBATCH --job-name=DMRG
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/DMRG-%A_%a.out
#SBATCH --error=logs/DMRG-%A_%a.err
#SBATCH --requeue

exec=${1}
inputfile=${2}
mkdir -p logs
ulimit -c unlimited
echo "Running job $SLURM_ARRAY_JOB_ID, iter $SLURM_ARRAY_TASK_ID at $HOSTNAME with inputfile $inputfile"
$exec $inputfile $SLURM_ARRAY_TASK_ID


