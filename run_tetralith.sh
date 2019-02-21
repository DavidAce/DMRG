#! /bin/bash

#SBATCH --exclusive
#SBATCH --job-name               DMRG
#SBATCH --time                   0-1:30:00
#SBATCH --time-min               0-0:30:00
#SBATCH --error                  log.error
#SBATCH --output                 log.out
#SBATCH --mem-per-cpu            2G
#SBATCH --cpus-per-task          1
#SBATCH --ntasks                 1
#SBATCH --ntasks-per-core        1
#SBATCH --kill-on-invalid-dep    yes


echo "Running on tetralith"
echo $@
exec=$1
file=$2
module load buildenv-gcc/2018a-eb
module load GCCcore/7.3.0
module load CMake/3.12.1
module load zlib/1.2.8
export CC=gcc
export CXX=g++
echo $LD_LIBRARY_PATH
LD_LIBRARY_PATH=$(gcc -print-file-name=libstdc++.so)
export LD_LIBRARY_PATH
export OMP_NUM_THREADS=1


$exec $file

