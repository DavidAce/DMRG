#!/bin/bash

echo "Running on tetralith"
echo $@
exec=$1
infile=$2
logfile=$3
errfile=$4

#SBATCH --exclusive
#SBATCH --job-name               DMRG
#SBATCH --time                   0-1:30:00
#SBATCH --time-min               0-0:30:00
#SBATCH --error                  $logfile
#SBATCH --output                 $errfile
#SBATCH --mem-per-cpu            2G
#SBATCH --cpus-per-task          1
#SBATCH --ntasks                 1
#SBATCH --ntasks-per-core        1
#SBATCH --kill-on-invalid-dep    yes

$exec $infile
