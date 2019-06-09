#! /bin/bash
# Run this file on tetralith with sbatch submit_slurm.sh
# This file will submit one job per file in ../input/

module load buildenv-gcc/2018a-eb
module load zlib/1.2.8
module load gimkl
module load GCCcore/7.3.0
module load clang/6.0.1
module load CMake/3.12.1
module load parallel/20181122-nsc1
source activate dmrg
#export CC=gcc
#export CXX=g++
export CC=clang
export CXX=clang++

export OMP_NUM_THREADS=1


exec=../build/Release/DMRG++
bunchfiles=$(find -L bunch -type f -name '*.txt')

for bunchfile in $bunchfiles; do
    [ -e "$bunchfile" ] || continue
    bunch_id=$(basename $bunchfile .txt)
    logname=logs/log_$bunch_id.out
    errname=logs/log_$bunch_id.err
    sbatch --output=$logname --error=$errname run_tetralith_bunch_parallel.sh $exec $bunchfile
done
