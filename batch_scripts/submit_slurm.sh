#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-b <build type>                     : Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-e                                  : Enable --exclusive mode. (default off)
-g                                  : Enable GNU Parallel. (default off)
-j <job name>                       : Job name. (default=DMRG)
-k <step size>                      : Step size of job arrays,for use with GNU Parallel (default = 32)
-m <memory (MB)>                    : Reserved amount of ram for each task in MB. (default = 4000)
-n <num sims>                       : Number of simulations per input file (default 10)
-o <other>                          : Other options passed to sbatch
-p <partition>                      : Partition name (default = all)
-r <requeue>                        : Enable --requeue, for requeuing in case of failure (default OFF)
-s <max simultaneous tasks>         : Simultaneous tasks for each job array (default = 32)
-t <time>                           : Time for each run (default = 1:00:00, i.e. 1 hour)
EOF
  exit 1
}

build=Release
nsims=10
maxtasks=32
jobname=DMRG
stepsize=32
mem=4000
time=--time=0-1:00:00
while getopts hb:egj:k:m:n:o:p:rs:t: o; do
    case $o in
        (h) usage ;;
        (b) build=$OPTARG;;
        (e) exclusive=--exclusive;;
        (g) gnuparallel=true;;
        (j) jobname=$OPTARG;;
        (k) stepsize=$OPTARG;;
        (m) mem=$OPTARG;;
        (n) nsims=$OPTARG;;
        (o) other=$OPTARG;;
        (p) partition=--partition=$OPTARG;;
        (r) requeue=--requeue;;
        (s) maxtasks=$OPTARG;;
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

while [ $((stepsize % nsims)) != 0  ]; do
    stepsize=$((stepsize-1))
    echo "adjusted stepsize = $stepsize"
done

inputfiles=$(find -L input -type f -name '*.cfg')
filecount=0
for inputfile in $inputfiles; do
    [ -e "$inputfile" ] || continue
    seedmin=$((filecount*nsims))
    seedmax=$((filecount*nsims+nsims-1))
    echo "Submitting jobs=[$seedmin - $seedmax]"

    if [ "$gnuparallel" = true ]; then
        module load parallel
        module load parallel/20181122-nsc1
        stepsize=$(( stepsize < nsims ? stepsize : nsims ))
        sbatch $partition $requeue $exclusive $time $other \
            --mem-per-cpu=$mem \
            --array=$seedmin-$seedmax:$stepsize --job-name=$jobname \
            run_jobarray_parallel.sh $exec $inputfile
    else
        sbatch $partition $requeue $exclusive $time $other \
            --mem-per-cpu=$mem \
            --array=$seedmin-$seedmax%$maxtasks --job-name=$jobname \
            run_jobarray.sh $exec $inputfile
    fi
    filecount=$((filecount+1))


done