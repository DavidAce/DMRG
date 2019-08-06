#!/bin/bash

read -p 'Number of realizations: ' -i 20 -e nsims
echo "Running $nsims realizations for each inputfile"

module load CMake
module load imkl
module load OpenBLAS
module load XZ/5.2.4-GCCcore-8.2.0
module load arpack-ng
module load ARPACK++
module load HDF5/1.10.5-GCCcore-8.2.0
module load Eigen
module load gflags
module load glog
module load CMake
module list

export OMP_NUM_THREADS=1


exec=../build/Release/DMRG++
if test -f "$exec"; then
    echo "Found exec: $exec"
else
    echo "exec does not exist: $exec"
    exit 1
fi

inputfiles=$(find -L input -type f -name '*.cfg')
count=0
nmax=$((nmin+nsims))
for inputfile in $inputfiles; do
    [ -e "$inputfile" ] || continue
    nmin=$((count*nsims))
    nmax=$((count*nsims+nsims-1))
    echo "Submitting array=[$nmin - $nmax]"
    sbatch  --array=$nmin-$nmax run_theophys_bunch_parallel.sh $exec $inputfile
    count=$((count+1))
done