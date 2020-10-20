#!/bin/bash
# Run this file on a slurm cluster, ./submit_slurm.sh


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-option | --option ] [=]<argument>
-h | --help                        : Help. Shows this text.
-b | --build-type <str>            : Build type: [ Release | RelWithDebInfo | Debug | Profile ]  (default = Release)
-c | --cluster <"a1,a2,a3...">     : Comma-separated string of clusters. At theophys: kraken | draken (default: )
     --cpus-per-task <int>         : Numer of cpu cores required per task (e.g. openmp threads) (default: 1)
-d | --dry-run                     : Do not actually submit
-e | --exclusive                   : Use nodes in exclusive node
     --exec-name <str>             : Name of the executable to run (default = DMRG++)
-f | --config <path or file>       : File or path to files containing simulation config files (suffixed .cfg) (default: input/ )
-J | --job-name <str>              : Slurm job name. (default = DMRG)
-m | --mem-per-cpu <int[suffix]>   : Memory per cpu core, e.g. 2000, 2000M, 2G (default: 2000)
-N | --sims-per-cfg <int>          : Number of simulations per config file. Can be split up into chunks with -n (default: 10)
-n | --sims-per-sbatch <int>       : Number of simulations per invocation of sbatch, or job-array size. (Splits -N into chunks of <num>)  (default = 10)
   | --ntasks <int>                : Number of tasks per simulation (e.g. mpi threads) (default: 1)
-o | --other <...>                 : Other options passed to sbatch
-O | --open-mode <append|truncate> : Open mode for logs (default: append)
-p | --partition <str>             : Partition name (default: dedicated)
-q | --qos <str>                   : Select Quality Of Service (default = )
-r | --requeue                     : Enable --requeue, for requeuing in case of failure (default OFF)
-s | --seed <path or file>         : File or path to files containing a list of seeds (suffixed .seed). Basenames should match the corresponding input files in -f (default: )
     --start-seed    <int>         : Starting seed (default: 0)
     --shuffle                     : Shuffle all seeeds and config files
-t | --time <time spec>            : Time for each run (default = 1:00:00, i.e. 1 hour)
-v | --verbose                     : Sets verbosity for sbatch

EOF
  exit 1
}

# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"   -o hb:c:def:J:m:N:n:O:p:q:rs:t:v \
                --long "\
                help\
                build-type:\
                cluster:\
                cpus-per-task:\
                dry-run\
                exclusive\
                exec-name:\
                config:\
                job-name:\
                mem-per-cpu:\
                sims-per-cfg:\
                sims-per-sbatch:\
                ntasks:\
                other:\
                open-mode:\
                partition:\
                qos:\
                requeue:\
                seed:\
                start-seed:\
                shuffle\
                time:\
                verbose\
                "  -- "$@")


#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ]; then exit 1 ; fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

execname=DMRG++
build_type=Release
jobname="--job-name=DMRG"
configpath=input/
startseed=0
time="--time=0-1:00:00"
simspercfg=10
simspersbatch=10
openmode="--open-mode=append"
cpuspertask="--cpus-per-task=1"
ntasks="--ntasks=1"

# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
#$1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
echo "Enabled options:"
while true;
do
  case "$1" in
    -h|--help)                      usage                                                                            ; shift   ;;
    -b|--build-type)                build_type=$2                     ; echo " * Build type               : $2"      ; shift 2 ;;
    -c|--cluster)                   cluster="--cluster=$2"            ; echo " * Cluster                  : $2"      ; shift 2 ;;
       --cpus-per-task)             cpuspertask="--cpus-per-task=$2"  ; echo " * Cpus per task            : $2"      ; shift 2 ;;
    -d|--dry-run)                   dryrun="ON"                       ; echo " * Dry run                  : ON"      ; shift   ;;
    -e|--exclusive)                 exclusive=--exclusive             ; echo " * Exclusive                : ON"      ; shift   ;;
       --exec-name)                 execname=$2                       ; echo " * Executable name          : $2"      ; shift 2 ;;
    -f|--config)                    configpath=$2                     ; echo " * Configs                  : $2"      ; shift 2 ;;
    -J|--job-name)                  jobname="--job-name=$2"           ; echo " * Job name                 : $2"      ; shift 2 ;;
    -m|--mem-per-cpu)               mempercpu="--mem-per-cpu=$2"      ; echo " * Mem per cpu              : $2"      ; shift 2 ;;
    -N|--sims-per-cfg)              simspercfg=$2                     ; echo " * Sims per cfg             : $2"      ; shift 2 ;;
    -n|--sims-per-sbatch)           simspersbatch=$2                  ; echo " * Sims per sbatch          : $2"      ; shift 2 ;;
       --ntasks)                    ntasks="--ntasks=$2"              ; echo " * ntasks                   : $2"      ; shift 2 ;;
       --other)                     other=$2                          ; echo " * Other                    : $2"      ; shift 2 ;;
    -O|--open-mode)                 openmode="--open-mode=$2"         ; echo " * Open mode                : $2"      ; shift 2 ;;
    -p|--partition)                 partition="--partition=$2"        ; echo " * Partition                : $2"      ; shift 2 ;;
    -q|--qos)                       qos="--qos=$2"                    ; echo " * Quality of service (QOS) : $2"      ; shift 2 ;;
    -r|--requeue)                   requeue="--requeue"               ; echo " * Requeue                  : ON"      ; shift   ;;
    -s|--seed)                      seedpath=$2                       ; echo " * Seed path                : $2"      ; shift   ;;
       --start-seed)                startseed=$2                      ; echo " * Start seed               : $2"      ; shift 2 ;;
       --shuffle)                   shuffle="ON"                      ; echo " * Shuffle                  : ON"      ; shift   ;;
    -t|--time)                      time="--time=0-$2"                ; echo " * Time                     : $2"      ; shift 2 ;;
    -v|--verbose)                   verbosity="-v"                    ; echo " * Verbosity                : ON"      ; shift   ;;
    --) shift; break;;
  esac
done

if [ $OPTIND -eq 0 ] ; then
    echo "No flags were passed"; usage ;exit 1;
fi
shift "$((OPTIND - 1))"

if [ "$sims-per-sbatch" -gt "$sims-per-cfg" ]; then
    echo "Cannot have sims-per-sbatch ($sims-per-sbatch) greater than sims-per-cfg ($sims-per-cfg)"
fi


exec=../build/$build_type/$execname
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
  rm -rf seeds
  mkdir -p seeds
  if [ -n "$shuffle" ]; then
      # Generate a master list with 2 columns: configs and seeds
      echo "Generating master file with config paths and seeds..."
      superfile=seeds/superfile.txt
      touch $superfile
      seedcounter=$((startseed + 0))
      for configfile in $configfiles; do
          for sim in $(seq $simspercfg); do
              echo "$configfile $seedcounter" >> $superfile
              seedcounter=$((seedcounter+1))
          done
      done
      echo "Generating master file with config paths and seeds... OK"

  else
      # Generate one list per cfg, each with 2 columns: config and seeds
      echo "Generating lists of config paths and seeds..."
      seedcounter=$((startseed + 0))
      for configfile in $configfiles; do
        configbase=$(basename $configfile .cfg)
        seedfile=seeds/$configbase.seed
        touch $seedfile
        for sim in $(seq $simspercfg); do
              echo "$configfile $seedcounter" >> $seedfile
              seedcounter=$((seedcounter+1))
          done
      done
      echo "Generating lists of config paths and seeds... OK"
  fi
fi


# Take seeds from seedfiles, distribute them to config files, randomize and then split into .seed files
if [ -n "$configfiles" ] && [ -n  "$seedfiles" ] && [ -n "$shuffle" ]; then
    rm -rf seeds
    mkdir -p seeds
    # Generate a master list with 2 columns: configs and seeds
    echo "Generating master file with config paths and given seeds..."
    superfile=seeds/superfile.txt
    touch $superfile
    for seedfile in $seedfiles; do
        seedbase=$(basename $seedfile .seed)
        configfile=$(find -L $configpath -type f -name $seedbase.cfg |  sort -g )
        if [ -z "$configfile" ]; then echo "Could not find a matching config file for seed file [ $seedfile ]. Searched with: find -L $configpath -type f -name $seedbase.cfg | Found: $match" ; exit 1; fi
        for seed in $(cat $seedfile); do
            echo "$configfile $seed" >> $superfile
        done
    done
    echo "Generating master file with config paths and given seeds... OK"
fi



if [ -e "$superfile" ] ; then
    echo "Shuffling and splitting master file into simulation files..."
    shuffledfile=seeds/randomsims.txt
    cat $superfile | shuf --output $shuffledfile
    split --lines=$simspersbatch --additional-suffix=.job -d --suffix-length=3 $shuffledfile seeds/part-
    rm $superfile $shuffledfile
    echo "Shuffling and splitting master file into simulation files... OK"
else
    echo "Splitting seed files into simulation files..."
    seedfiles=$(find -L seeds -type f -name '*.seed' |  sort -g )
    for seedfile in $seedfiles; do
      seedbase=$(basename $seedfile .seed)
      split --lines=$simspersbatch --additional-suffix=.job -d --suffix-length=3 $seedfile seeds/$seedbase-
    done


    echo "Splitting seed files into simulation files... OK"
fi




# From this point on we are guaranteed to have a set of files in seeds/part_###.job or seeds/<configbase>-###.job
# Each .job file contains 2 columns with the corresponding config files and seeds


jobfiles=$(find -L seeds -type f -name '*.job' |  sort -g )
if [ -z "$jobfiles" ] ; then
    echo "No simulation files found!"
    exit 1
fi

for jobfile in $jobfiles; do
  if [ -n "$dryrun" ]; then
cat << EOF >&2
sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask \
--array=1-$simspersbatch \
run_jobarray.sh -e $exec -f $jobfile
EOF
  bash run_jobarray.sh -e $exec -f $jobfile -d
  else
    sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask \
      --array=1-$simspersbatch \
      run_jobarray.sh -e $exec -f $jobfile
  fi
done

