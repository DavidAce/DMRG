#! /bin/bash
# Run this file on tetralith with sbatch submit_slurm.sh
# This file will submit one job per file in ../input/

read -p 'Number of realizations [100]: ' -e -i 100 nsims
echo 'Running $nsims realizations for each inputfile'


module load buildenv-gcc/2018a-eb
module load zlib/1.2.8
module load CMake/3.12.1
source activate dmrg
module load gimkl
module load GCCcore/7.3.0
module load clang/6.0.1
module load CMake/3.12.1





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


#
#exec=../build/Release/DMRG++
#bunchfiles=$(find -L bunch -type f -name '*.txt')
#
#for bunchfile in $bunchfiles; do
#    [ -e "$bunchfile" ] || continue
#    bunch_id=$(basename $bunchfile .txt)
#    logname=logs/log_$bunch_id.out
#    errname=logs/log_$bunch_id.err
#    sbatch --output=$logname --error=$errname run_tetralith_bunch_parallel.sh $exec $bunchfile
#done
