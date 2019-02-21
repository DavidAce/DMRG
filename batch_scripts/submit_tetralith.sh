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


base=mbl_
exec=../build/Release/DMRG++

for inputfile in ../input/L_*/$(NAME)_*.cfg; do
    [ -e "$inputfile" ] || continue
    base_id=$(basename $inputfile .cfg)
    logname=logs/log_$base_id.out
    errname=logs/log_$base_id.err
    sbatch run_tetralith.sh $exec $inputfile $logname $errname
done
