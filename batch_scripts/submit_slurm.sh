#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                : $PROGNAME [-options] with the following options:
-h                   : Help. Shows this text.
-b <build type>      : Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-e                   : Enable --exclusive mode. (default off)
-j <job name>        : Job name. (default=DMRG)
-m <memory (MB)>     : Reserved amount of ram for each task in MB. (default = 4000)
-n <num sims>        : Number of simulations (default 10)
-o <other>           : Other options passed to sbatch
-p <partition>       : Partition name (default = all)
-r <requeue>         : Enable --requeue, for requeuing in case of failure (default OFF)
-s <num simult.>     : Max number of simultaneous tasks to run for each job array (default = 32)
-t <time>            : Time for each run (default = 1:00:00, i.e. 1 hour)
EOF
  exit 1
}

build=Release
nsims=10
nsimu=%32
jobname=DMRG
mem=4000
time=--time=0-1:00:00
while getopts hb:ej:m:n:o:p:rt: o; do
    case $o in
        (h) usage ;;
        (b) build=$OPTARG;;
        (e) exclusive=--exclusive;;
        (j) jobname=$OPTARG;;
        (m) mem=$OPTARG;;
        (n) nsims=$OPTARG;;
        (o) other=$OPTARG;;
        (p) partition=--partition=$OPTARG;;
        (r) requeue=--requeue;;
        (s) nsimu=%$OPTARG;;
        (t) time=--time=0-$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


echo "Running $nsims realizations for each inputfile"

export OMP_NUM_THREADS=1


exec=../build/$build/DMRG++
if test -f "$exec"; then
    echo "Found exec: $exec"
else
    echo "exec does not exist: $exec"
    exit 1
fi

inputfiles=$(find -L input -type f -name '*.cfg')
count=0
for inputfile in $inputfiles; do
    [ -e "$inputfile" ] || continue
    nmin=$((count*nsims))
    nmax=$((count*nsims+nsims-1))
    #echo "Submitting array=[$nmin - $nmax]"
    sbatch $partition $requeue $exclusive $time $other \
            --hint=compute_bound \
            --array=$nmin-$nmax$nsimu --job-name=$jobname \
            run_jobarray.sh $exec $inputfile
    count=$((count+1))
done