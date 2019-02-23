#! /bin/bash
# Run this file on tetralith with sbatch submit_slurm.sh
# This file will submit one job per file in ../input/


module load buildenv-gcc/2018a-eb
module load GCCcore/7.3.0
export CC=gcc
export CXX=g++
echo $LD_LIBRARY_PATH
LD_LIBRARY_PATH=$(gcc -print-file-name=libstdc++.so)
export LD_LIBRARY_PATH
export OMP_NUM_THREADS=1


exec=../build/Release/DMRG++
bunchfiles=$(find -L bunch -type f -name '*.txt')

for bunchfile in $bunchfiles; do
    [ -e "$bunchfile" ] || continue
    bunch_id=$(basename $bunchfile .txt)
    logname=logs/log_$bunch_id.out
    errname=logs/log_$bunch_id.err
    sbatch --output=$logname --error=$errname run_tetralith_bunch.sh $exec $bunchfile
done
