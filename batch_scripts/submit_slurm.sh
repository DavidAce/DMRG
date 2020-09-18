#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-a <app name>                       : Name of the executable to run (default = DMRG++)
-b <build type>                     : Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-c <"cluster_list">                 : String with comma-separated list of clusters. At theophys: kraken | draken (default = )
-e                                  : Enable --exclusive mode. (default = off)
-f <config file/path>               : File or path to files containing simulation config files (suffixed .cfg) (default = input/ )
-g <seeds  file/path>               : File or path to files containing a list of seeds (suffixed .seed). Basenames should match the corresponding input files in -f (default = )
-j                                  : Enable job-array mode with a certain array size.(default = "OFF")
-J <job name>                       : Job name. (default = DMRG)
-m <memory (MB)>                    : Reserved amount of ram for each task in MB.
-n <ntasks (sims in parallel)>      : Number of simulations that are run in parallel by GNU Parallel (set to 1 for job arrays) (default = 32)
-N <num sims per cfg>               : Number of simulations per config file (default = 10)
-o <other>                          : Other options passed to sbatch
-p <partition>                      : Partition name (default = dedicated)
-q <qos>                            : Select Quality Of Service (default = )
-r <requeue>                        : Enable --requeue, for requeuing in case of failure (default OFF)
-s <sims per sbatch>                : Number of tasks/simulations (or job-array size) per invocation of sbatch  (default = 10)
-S <start seed>                     : Starting seed, if you don't want to start from 0.
-t <time>                           : Time for each run (default = 1:00:00, i.e. 1 hour)
-v <vv...>                          : Sets verbosity for sbatch
EOF
  exit 1
}

execname=DMRG++
build=Release
jobname=DMRG
configpath=input/
startseed=0
time=--time=0-1:00:00
ntasks_parallel=32
simspercfg=10
simspersbatch=10

while getopts ha:b:c:def:g:jJ:m:n:N:o:p:rs:S:t:v: o; do
    case $o in
        (h) usage ;;
        (a) execname=$OPTARG;;
        (b) build=$OPTARG;;
        (c) cluster="--cluster=$OPTARG";;
        (d) dryrun=true;;
        (e) exclusive=--exclusive;;
        (f) configpath=$OPTARG;;
        (g) seedpath=$OPTARG;;
        (j) jobarray="ON";;
        (J) jobname=$OPTARG;;
        (m) mempercpu="--mem-per-cpu=$OPTARG";;
        (n) ntasks_parallel=$OPTARG;;
        (N) simspercfg=$OPTARG;;
        (o) other=$OPTARG;;
        (p) partition="--partition=$OPTARG";;
        (q) qos="--qos=$OPTARG";;
        (r) requeue=--requeue;;
        (s) simspersbatch=$OPTARG;;
        (S) startseed=$OPTARG;;
        (t) time="--time=0-$OPTARG";;
        (v) verbosity="-v$OPTARG";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


export OMP_NUM_THREADS=1


# Handle incompatible options
#if [ -n "$seedpath" ] ; then
#    if [ "$startseed" -gt 0 ] ; then echo "Options  -g and -S are incompatible" ; exit 1; fi
#    # We only need simspernode to continue
#    if   [ -z "$simspernode" ]; then echo "Given option -g, you must also give -s" ;exit 1; fi
#    if   [ -n "$simspercfg" ];  then echo "Given option -g, you must also give -s, but not -n" ;exit 1; fi
#    if   [ -n "$nodespercfg" ]; then echo "Given option -g, you must also give -s, but not -N" ;exit 1; fi
#    echo "Sims per node = $simspernode"
#
#else
#    if [ -z "$simspercfg" ]  && [ -z "$simspernode" ];  then echo "At least two of -N -n and -s must be given." ;exit 1; fi
#    if [ -z "$simspercfg" ]  && [ -z "$nodespercfg" ];  then echo "At least two of -N -n and -s must be given." ;exit 1; fi
#    if [ -z "$simspernode" ] && [ -z "$nodespercfg" ]; then echo "At least two of -N -n and -s must be given." ;exit 1; fi
#
#    if [ -n "$nodespercfg" ] && [ "$nodespercfg" -gt 1 ]; then
#        # Make sure simspercfg and simspernode haven't both been given
#        if [ -n "$simspercfg" ] && [ -n "$simspernode" ]; then echo "Options -N -n and -s are incompatible together. Choose two."  ; exit 1; fi
#        if [ -n "$simspercfg" ] && [ -z "$simspernode" ]; then echo "Computing simspernode" ; simspernode=$((simspercfg / nodespercfg)); fi
#        if [ -z "$simspercfg" ] && [ -n "$simspernode" ]; then echo "Computing simspercfg" ; simspercfg=$((simspernode * nodespercfg)); fi
#    else
#        # Make sure both simspercfg and simspernode is given,
#        if [ -n "$simspercfg" ] && [ -n "$simspernode" ]; then echo "Computing nodespercfg"; nodespercfg=$((simspercfg / simspernode)) ; fi
#    fi
#
#    echo "Sims per conf = $simspercfg"
#    echo "Sims per node = $simspernode"
#    echo "Node per conf = $nodespercfg"
#fi





#if [ "$simspernode" -lt 32 ]; then
#    read -p "Running with fewer than 32 simulations per node. Are you sure? [y/n] " -n 1 -r
#    echo    # (optional) move to a new line
#    if [[ $REPLY =~ ^[Yy]$ ]]; then
#        echo "Going ahead with $simspernode simulations per node"
#    else
#        echo "Abort."
#        exit 0
#    fi
#fi




exec=../build/$build/$execname
if [ -f "$exec" ]; then
    echo "Found executable: $exec"
else
    echo "Executable does not exist: $exec"
    exit 1
fi

# Load the required parallel module
if [[ "$HOSTNAME" == *"tetralith"* ]];then
    module try-load parallel/20181122-nsc1
else
    module try-load parallel
fi


# The point from now on will be to generate one massive list with 2 columns: configs and seeds
# After generating that list, we shuffle it, and then split it into chunks

# Parse config path
echo "Finding config path..."
if [ -d "$configpath" ];then
   configfiles=$(find -L $configpath -type f -name '*.cfg' |  sort -g )
   if [ -z "$configfiles" ]; then
        echo "No config files found"
        exit 1
    fi
elif [ -f "$configpath" ]; then
    configfiles=$configpath
    if [ ! -e "$configfiles" ]; then
        echo "Config file does not exist: $configfiles"
        exit 1
    fi
else
    echo "Config file path cannot be parsed: $configpath "
    exit 1
fi
echo "Finding config path... OK"


# Parse seed path
if [ -n "$seedpath" ] ; then
    if [ -d "$seedpath" ];then
        echo "Finding seed files in given path..."
        seedfiles=$(find -L $seedpath -type f -name '*.seed' |  sort -g )
        echo "Finding seed files in given path... OK"
    elif [ -f "$seedpath" ]; then
        echo "Found seed file in given path: $seedpath"
        seedfiles=$seedpath
        seedpath=$(dirname $seedpath)
    elif [ ! -e "$seedpath" ]; then
        echo "Seed path does not exist: $seedpath "
        exit 1
    else
        echo "Input file path cannot be parsed: $seedpath "
        exit 1
    fi
    echo "Pairing seed files with config files..."
    for seedfile in $seedfiles; do
        echo "Matching seed file: $seedfile"
        seedbase=$(basename $seedfile .seed)
        match=$(find -L $configpath -type f -name $seedbase.cfg |  sort -g )
        num=$(echo $match | wc -w)
        if [ -z "$match" ]  ; then echo "Could not find a config file matching the given seed file [ $seedfile ]. Searched with: find -L $configpath -type f -name '$seedbase.cfg' | Found: $match" ; exit 1; fi
        if [ "$num" -gt 1 ] ; then echo "Too many config files correspond to seed file [ $seedfile ]. Found $num files: $match"; exit 1; fi
    done
    echo "Pairing seed files with config files... OK"
fi



# Generate seeds, distribute them to config files, randomize and then split into .seed files
if [ -n "$configfiles" ] && [ -z  "$seedfiles" ] ; then
    # Generate a master list with 2 columns: configs and seeds
    echo "Generating master file with config paths and seeds..."
    masterfile=seeds/masterfile.txt
    rm -rf seeds
    mkdir -p seeds
    touch $masterfile
    seedcounter=$((startseed + 0))
    for configfile in $configfiles; do
        for sim in $(seq $simspercfg); do
            echo "$configfile $seedcounter" >> $masterfile
            seedcounter=$((seedcounter+1))
        done
    done
    echo "Generating master file with config paths and seeds... OK"
fi




# Take seeds from seedfiles, distribute them to config files, randomize and then split into .seed files
if [ -n "$configfiles" ] && [ -n  "$seedfiles" ] ; then
    # Generate a master list with 2 columns: configs and seeds
    echo "Generating master file with config paths and given seeds..."
    masterfile=seeds/masterfile.txt
    rm -rf seeds
    mkdir -p seeds
    touch $masterfile
    for seedfile in $seedfiles; do
        seedbase=$(basename $seedfile .seed)
        configfile=$(find -L $configpath -type f -name $seedbase.cfg |  sort -g )
        if [ -z "$configfile" ]; then echo "Could not find a matching config file for seed file [ $seedfile ]. Searched with: find -L $configpath -type f -name $seedbase.cfg | Found: $match" ; exit 1; fi
        for seed in $(cat $seedfile); do
            echo "$configfile $seed" >> $masterfile
        done
    done
    echo "Generating master file with config paths and given seeds... OK"
fi



if [ -e "$masterfile" ] ; then
    echo "Shuffling and splitting master file into simulation files..."
    shuffledfile=seeds/randomsims.txt
    cat $masterfile | shuf --output $shuffledfile
    split --lines=$simspersbatch --additional-suffix=.sim -d --suffix-length=3 $shuffledfile seeds/part_
    rm $masterfile $shuffledfile
    echo "Shuffling and splitting master file into simulation files... OK"
else
    echo "Master file could not be found: $masterfile"
fi




# From this point on we are guaranteed to have a set of files in seeds/part_###.sim
# Each .sim file contains 2 randomized columns with config files and corresponding seeds


simfiles=$(find -L seeds -type f -name '*.sim' |  sort -g )
if [ -z "$simfiles" ] ; then
    echo "No simulation files found!"
    exit 1
fi

for simfile in $simfiles; do
    if [ -n "$dryrun" ] ; then
        echo "sbatch $cluster $partition $mempercpu $requeue $exclusive $time $other $verbosity --job-name=$jobname --ntasks=$ntasks_parallel run_parallel.sh -e $exec -f $simfile"
    elif [ -n "$jobarray" ]; then
        sbatch $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $verbosity \
            --job-name=$jobname \
            --ntasks=$ntasks_parallel \
            --array=1-$simspersbatch \
            run_jobarray.sh -e $exec -f $simfile
    else
        sbatch $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $verbosity \
            --job-name=$jobname \
            --ntasks=$ntasks_parallel \
            run_parallel.sh -e $exec -f $simfile
    fi
done

