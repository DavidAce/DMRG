#!/bin/bash


#SBATCH --job-name               =DMRG
#SBATCH --time                   =0-1:30:00
#SBATCH --time-min               =0-0:30:00
#SBATCH --mem-per-cpu            =2000
#SBATCH --cpus-per-task          =1
#SBATCH --ntasks                 =1
#SBATCH --ntasks-per-core        =1
#SBATCH --kill-on-invalid-dep    =yes

echo "Running on tetralith"
echo $@
exec=$1
infile=$2
$exec $infile
