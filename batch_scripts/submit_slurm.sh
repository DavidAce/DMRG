#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage            : $PROGNAME  [-h ] [-m <mode>] [-n <num sims>] [-p <partition>] [-t <time>]
-h               : Help. Shows this text.
-m <mode>        : Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-n <num sims>    : Number of simulations (default 10)
-p <partition>   : Partition name (default = )
-t <time>        : Time for each run (default = 4:00:00)
EOF
  exit 1
}

mode=Release
nsims=10
partition=all
time=0-4:00:00
while getopts hm:n:p:t: o; do
    case $o in
        (h) usage ;;
        (m) mode=$OPTARG;;
        (n) nsims=$OPTARG;;
        (p) partition=$OPTARG;;
        (t) time=0-$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


echo "Running $nsims realizations for each inputfile"

export OMP_NUM_THREADS=1


exec=../build/$mode/DMRG++
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
    sbatch --partition=$partition --time=$time --array=$nmin-$nmax run_jobarray.sh $exec $inputfile
    count=$((count+1))
done