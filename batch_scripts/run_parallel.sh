#!/bin/bash

#SBATCH --exclusive
#SBATCH --job-name=DMRG
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --kill-on-invalid-dep=yes


exec=${1}
bunchfile=${2}
bunchbase=$(basename $bunchfile .txt)
outdir=logs/$bunchbase
mkdir -p $outdir

cat $bunchfile | parallel  $exec {} ">" $outdir/{/.}.out
