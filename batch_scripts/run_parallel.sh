#!/bin/bash

#SBATCH --nodes=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/DMRG-%A.out
#SBATCH --error=logs/DMRG-%A.err




PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-e <executable>                     : Path to executable (default = "")
-f <inputfile>                      : Path to inputfile with settings (default = "")
-i <input seed file>                : Path to file containing a list of seeds. Incompatible with -l and -u. (default = "")
-l <lower seed>                     : Lower bound of seed range (default = )
-u <upper seed>                     : Upper bound of seed range (default = )
-o <output directory>               : Output directory for logs (default = "")
EOF
  exit 1
}

exec=""
inputfile=""
seedfile=""
nmin=
nmax=

while getopts he:f:i:l:u:o: o; do
    case $o in
        (h) usage ;;
        (e) exec=$OPTARG;;
        (f) inputfile=$OPTARG;;
        (i) seedfile=$OPTARG;;
        (l) nmin=$OPTARG;;
        (u) nmax=$OPTARG;;
        (o) outdir=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done

if [ -n "$nmin" ]; then
    if [ -z "$nmax" ]; then
        echo "Both bounds are needed, -l and -u."
    elif [ -z "$outdir" ] ; then
        outdir=logs/seed_$nmin-$nmax
    fi
fi

if [ -n "$seedfile" ] ; then
    if [ -n "$nmin" ] || [ -n "$nmax" ]; then
        echo "Can't both have a range [nmin nmax] and a seedfile"
        exit 1
    elif [ -z "$outdir" ] ; then
        outdir=logs/seed_from_file
    fi
fi




inputbase=$(basename $inputfile .cfg)
mkdir -p $outdir

echo "CPUS ON  NODE  : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE  : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK  : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU    : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE   : $SLURM_MEM_PER_NODE"

if [ -n "$seedfile" ]; then
    echo "Running job $SLURM_JOB_ID, with seeds from file $seedfile at $HOSTNAME with inputfile $inputfile"
    cat $seedfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log "$exec -i $inputfile -r {} &> $outdir/${inputbase}_{}.out"
elif [ -n "$nmin" ] && [ -n "$nmax" ] ;then
    echo "Running job $SLURM_JOB_ID, seeds [$nmin - $nmax] at $HOSTNAME with inputfile $inputfile"
    parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log "$exec -i $inputfile -r {} &> $outdir/${inputbase}_{}.out" ::: $(seq $nmin $nmax)
fi


# For initial state enumeration
#realization=1811
#parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log "$exec -i $inputfile -r $realization -s {} &> $outdir/${inputbase}_${realization}_{}.out" ::: $(seq $nmin $nmax)






# parallel --memfree $SLURM_MEM_PER_CPU --joblog $outdir/$inputbase.log  $exec $inputfile {} '2>&1' ">" $outdir/$inputbase_{}.out ::: $(seq $nmin $nmax)

# Explanation

# $(seq $nmin $nmax) will print a sequence of seeds from nmin to nmax. Each number is piped with | to parallel,
# where the number is found in {1}.
# The lines are piped with "|" into GNU parallel, which takes a list as input argument. Each element of the list
# is expanded with {}.

# For each line, parallel will run exec with that line as input argument. The output of each execution is piped with ">"
# into file $outdir/$inputbase_{}.out   which will expand to something like  logs/jobarray_0/mbl_blabla_123.out

# There are many, many more arguments you can give to parallel, to control the number of cores, memory requirements per thread etc.