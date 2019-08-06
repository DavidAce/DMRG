#! /bin/bash
# Run this file on slurm cluster with sbatch submit_slurm.sh
# This file will submit one job per file in ../input/

read -p 'Number of realizations: ' -i 20 -e nsims
echo "Running $nsims realizations for each inputfile"



export OMP_NUM_THREADS=1
exec=../build/Release/DMRG++
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


