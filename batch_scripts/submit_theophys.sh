#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-c] [-h ] [-j <num_threads>] [-l] [-m <mode>] [-t <target>] [-a <march>]
-h               : Help. Shows this text.
-m <mode>        : Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-n <sims>        : Number of simulations (default 10)
EOF
  exit 1
}

mode="Release"
nsims=10
while getopts hm:n: o; do
    case $o in
        (h) usage ;;
        (m) mode=$OPTARG;;
        (n) nsims=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


echo "Running $nsims realizations for each inputfile"

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