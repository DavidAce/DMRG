#!/bin/bash

#SBATCH --exclusive
#SBATCH --job-name=DMRG256
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-36:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --kill-on-invalid-dep=yes


exec=${1}
bunchfile=${2}
bunchbase=$(basename $bunchfile .txt)
outdir=logs/$bunchbase
mkdir -p $outdir

cat $bunchfile | parallel --colsep ' ' $exec {1} {2} '2>&1' ">" $outdir/{1/.}_{2}.out


# Explanation

# cat $bunchfile will print all the lines in a bunchfile. Each line is in the bunchfile the path to an input file that
# the executable exec takes as input argument. The lines look like "../input/L_12/mbl_232_blabla.cfg"

# The lines are piped with "|" into GNU parallel, which takes a list as input argument. Each element of the list
# is expanded with {}, and you can extract from {} a basename without extension  using {/.} See more in the manual for
# parallel

# For each line, parallel will run exec with that line as input argument. The output of each execution is piped with ">"
# into file $outdir/{/.}.out   which will expand to something like  logs/bunch_0/mbl_232_blabla.out

# There are many, many more arguments you can give to parallel, to control the number of cores, memory requirements per thread etc.
