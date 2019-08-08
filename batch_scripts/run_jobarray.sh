#!/bin/bash
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/DMRG-%A_%a.out
#SBATCH --error=logs/DMRG-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1


exec=${1}
inputfile=${2}
mkdir -p logs
ulimit -c unlimited
echo "Running job $SLURM_ARRAY_JOB_ID, iter $SLURM_ARRAY_TASK_ID at $HOSTNAME with inputfile $inputfile"
echo "CPUS ON  NODE: $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE: $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK: $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU  : $SLURM_MEM_PER_CPU"
$exec $inputfile $SLURM_ARRAY_TASK_ID

