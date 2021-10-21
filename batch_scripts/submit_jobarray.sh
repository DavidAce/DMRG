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
     --ntasks-per-core             : Number of tasks on each core (default: 1)
-d | --dry-run                     : Do not actually submit
-e | --exclusive                   : Use nodes in exclusive node
     --execname <str>              : Name of the executable to run (default = DMRG++)
-f | --config <path or file>       : File or path to files containing simulation config files (suffixed .cfg) (default: input/ )
-j | --job-file <path or file>     : File containing config-seed pairs. Use for resuming failed runs
-J | --job-name <str>              : Slurm job name. (default = DMRG)
-m | --mem-per-cpu <int[suffix]>   : Memory per cpu core, e.g. 2000, 2000M, 2G (default: 2000)
-N | --sims-per-cfg <int>          : Number of simulations per config file. Can be split up into chunks with -n (default: 10)
-n | --sims-per-array <int>        : Number of simulations in each job-array (splits -N into -N/-n chunks)  (default: 1000)
   | --sims-per-task <int>         : Number of simulations per job-array task. This is equivalent to the step, or stride in the array (default: 1)
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
PARSED_OPTIONS=$(getopt -n "$0"   -o hb:c:def:j:J:m:N:n:o:O:p:q:r:s:t:v \
                --long "\
                help\
                build-type:\
                cluster:\
                cpus-per-task:\
                ntasks-per-core:\
                dry-run\
                exclusive\
                execname:\
                config:\
                job-file:\
                job-name:\
                mem-per-cpu:\
                sims-per-cfg:\
                sims-per-array:\
                sims-per-task:\
                ntasks:\
                other:\
                open-mode:\
                partition:\
                qos:\
                requeue\
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
simsperarray=1000
simspertask=1
openmode="--open-mode=append"
cpuspertask="--cpus-per-task=1"
ntaskspercore="--ntasks-per-core=1"
ntasks="--ntasks=1"

# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
#$1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
echo "Enabled options:"
while true;
do
  case "$1" in
    -h|--help)                      usage                                                                            ; shift   ;;
    -b|--build-type)                build_type=$2                         ; echo " * Build type               : $2"      ; shift 2 ;;
    -c|--cluster)                   cluster="--cluster=$2"                ; echo " * Cluster                  : $2"      ; shift 2 ;;
       --cpus-per-task)             cpuspertask="--cpus-per-task=$2"      ; echo " * Cpus per task            : $2"      ; shift 2 ;;
       --ntasks-per-core)           ntaskspercore="--ntasks-per-core=$2"  ; echo " * Tasks per core           : $2"      ; shift 2 ;;
    -d|--dry-run)                   dryrun="ON"                           ; echo " * Dry run                  : ON"      ; shift   ;;
    -e|--exclusive)                 exclusive=--exclusive                 ; echo " * Exclusive                : ON"      ; shift   ;;
       --execname)                  execname=$2                           ; echo " * Executable name          : $2"      ; shift 2 ;;
    -f|--config)                    configpath=$2                         ; echo " * Configs                  : $2"      ; shift 2 ;;
    -j|--job-file)                  jobfile=$2                            ; echo " * Job file                 : $2"      ; shift 2 ;;
    -J|--job-name)                  jobname="--job-name=$2"               ; echo " * Job name                 : $2"      ; shift 2 ;;
    -m|--mem-per-cpu)               mempercpu="--mem-per-cpu=$2"          ; echo " * Mem per cpu              : $2"      ; shift 2 ;;
    -N|--sims-per-cfg)              simspercfg=$2                         ; echo " * Sims per cfg             : $2"      ; shift 2 ;;
    -n|--sims-per-array)            simsperarray=$2                       ; echo " * Sims per array           : $2"      ; shift 2 ;;
       --sims-per-task)             simspertask=$2                        ; echo " * Sims per task            : $2"      ; shift 2 ;;
       --ntasks)                    ntasks="--ntasks=$2"                  ; echo " * ntasks                   : $2"      ; shift 2 ;;
     o|--other)                     other=$2                              ; echo " * Other                    : $2"      ; shift 2 ;;
    -O|--open-mode)                 openmode="--open-mode=$2"             ; echo " * Open mode                : $2"      ; shift 2 ;;
    -p|--partition)                 partition="--partition=$2"            ; echo " * Partition                : $2"      ; shift 2 ;;
    -q|--qos)                       qos="--qos=$2"                        ; echo " * Quality of service (QOS) : $2"      ; shift 2 ;;
    -r|--requeue)                   requeue="--requeue"                   ; echo " * Requeue                  : ON"      ; shift   ;;
    -s|--seed)                      seedpath=$2                           ; echo " * Seed path                : $2"      ; shift 2 ;;
       --start-seed)                startseed=$2                          ; echo " * Start seed               : $2"      ; shift 2 ;;
       --shuffle)                   shuffle="ON"                          ; echo " * Shuffle                  : ON"      ; shift   ;;
    -t|--time)                      time="--time=0-$2"                    ; echo " * Time                     : $2"      ; shift 2 ;;
    -v|--verbose)                   verbosity="-v"                        ; echo " * Verbosity                : ON"      ; shift   ;;
    --) shift; break;;
  esac
done

if [ $OPTIND -eq 0 ] ; then
    echo "No flags were passed"; usage ;exit 1;
fi
shift "$((OPTIND - 1))"


if [ "$simsperarray" -gt "$simspercfg" ]; then
    echo "Cannot have array-size ($simsperarray) > sims-per-cfg ($simspercfg)"
    exit 1
fi

# Deactivate conda so that we do not get conflicting dynamic libraries
if [ -n "$CONDA_PREFIX" ] ; then
    if [ -f "$CONDA_PREFIX_1/etc/profile.d/conda.sh" ]; then
        source $CONDA_PREFIX_1/etc/profile.d/conda.sh
    fi
    if [ -f "$CONDA_PREFIX/etc/profile.d/conda.sh" ]; then
        source $CONDA_PREFIX/etc/profile.d/conda.sh
    fi
    counter=0
    while [ -n "$CONDA_PREFIX" ] && [  $counter -lt 10 ]; do
        let counter=counter+1
        echo "Deactivating conda environment $CONDA_PREFIX"
        conda deactivate
    done
    conda info --envs
fi


exec=$(find ../build/$build_type ../build/$build_type/bin -maxdepth 1 -type f -executable -name "$execname" -print -quit)
if [ -f "$exec" ]; then
    echo "Found executable $execname: $exec"
else
    echo "Executable $execname does not exist: $exec"
    exit 1
fi

# Make sure the executable would find its dynamically loaded libraries
# This is copied from build.sh
if [[ "$HOSTNAME" == *"tetralith"* ]];then
    echo "Running on tetralith"
    module load foss/2020b
    export MKLROOT=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a/mkl
    export EBROOTIMKL=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a
elif [[ "$HOSTNAME" == *"raken"* ]];then
    module load imkl
fi

if ldd $exec | grep 'not found'; then
  echo "Some dynamic libraries were not found."
  echo "Perhaps a module needs to be loaded."
  exit 1
fi




# Load the required parallel module
if [[ "$HOSTNAME" == *"tetralith"* ]];then
    module try-load parallel/20181122-nsc1
else
    module try-load parallel
fi


# If we are using a job-file we can use it directly and then skip generating one
#      split --lines=$simsperarray --additional-suffix=.job -d --suffix-length=3 $seedfile seeds/$seedbase-

if [ -f "$jobfile" ]; then
  numseeds=$(cat $jobfile | wc -l)
  jobbase=$(basename $jobfile .job)
  jobdir=$(dirname $jobfile)
  if [ "$numseeds" -gt $simsperarray ] ; then
    split --lines=$simsperarray --additional-suffix=.job -d --suffix-length=3 $jobfile $jobdir/$jobbase-
    mv "$jobfile" "$(basename $jobfile).bak"
  fi

  jobfiles=$(find -L $jobdir -type f -name '*.job' |  sort -g )
  if [ -z "$jobfiles" ] ; then
      echo "No simulation files found!"
      exit 1
  fi

  if [ -z "$dryrun" ]; then
    git log -1 >> job_report.txt
  fi

  for file in $jobfiles; do
    numseeds=$(cat $file | wc -l)

    if [ -n "$dryrun" ]; then
cat << EOF >&2
sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
--array=1-$numseeds:$simspertask \
run_jobarray.sh -e $exec -f $file
EOF
    bash run_jobarray.sh -e $exec -f $file -d
    else
      echo "sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
        --array=1-$numseeds:$simspertask \
        run_jobarray.sh -e $exec -f $file" >> job_report.txt

      sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
        --array=1-$numseeds:$simspertask \
        run_jobarray.sh -e $exec -f $file
    fi
  done
  exit 0
fi


# The point from now on will be to generate one massive list with 2 columns: configs and seeds
# After generating that list we split it into chunks and write the chunks to .job files
# Optionally the big list is shuffled before the split

# Parse config path
if [ -d "$configpath" ];then
   echo "Finding config path..."
   configfiles=$(find -L $configpath -type f -name '*.cfg' |  sort -g )
   if [ -z "$configfiles" ]; then
        echo "No config files found"
        exit 1
    fi
elif [ -f "$configpath" ]; then
    echo "Finding config path..."
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


# Generate a new unique seed path
date=$(date '+%Y-%m-%dT%H:%M:%S')
seedpath="seeds-$date"
if [ -d $seedpath ]; then
  n=1
  seedpath="seeds-$date-$n"
  # Increment $n as long as a directory with that name exists
  while [[ -d "$seedpath" ]] ; do
      n=$(($n+1))
      seedpath="seeds-$date-$n"
  done
fi
mkdir -p $seedpath


# Generate seeds, distribute them to config files, optionally shuffle, and then split into .seed files
if [ -n "$configfiles" ] && [ -z  "$seedfiles" ] ; then
  echo "Generating list of config paths and seeds..."
  if [ -n "$shuffle" ]; then
      # Generate a list with 2 columns: configs and seeds
      # The superlist has length simspercfg * num config files
      # The superlist is shuffled then split back into lists of length simspercfg

      superfile=$seedpath/superfile.seed
      truncate -s 0 $superfile
      seedcounter=$((startseed + 0))
      for configfile in $configfiles; do
          for sim in $(seq $simspercfg); do
              echo "$configfile $seedcounter" >> $superfile
              seedcounter=$((seedcounter+1))
          done
      done
      echo "Shuffling seeds ..."
      shuffledfile=$seedpath/shuffledfile.txt
      truncate -s 0 $shuffledfile
      cat $superfile | shuf --output $shuffledfile
      echo "Shuffling seeds ... OK"
      echo "Splitting list ..."
      seedcounter=$((startseed + 0))
      for configfile in $configfiles; do
          endseed=$((seedcounter+simspercfg-1))
          seedfile="$seedpath/part.[$seedcounter.$endseed].seed"
          touch $seedfile
          # Move simspercfg lines from shuffledfile --> seedfile
          head -$simspercfg $shuffledfile > $seedfile && sed -i '1,$simspercfgd' $shuffledfile
          seedcounter=$((seedcounter+simspercfg))
      done
      rm $superfile $shuffledfile
      echo "Splitting list ... OK"
  else
      # Generate a list with 2 columns: configs and seeds.
      # Each list has length simspercfg
      seedcounter=$((startseed + 0))
      for configfile in $configfiles; do
          configbase=$(basename $configfile .cfg)
          endseed=$((seedcounter+simspercfg-1))
          seedfile="$seedpath/$configbase.[$seedcounter-$endseed].seed"
          tempfile="/dev/shm/$configbase.[$seedcounter-$endseed].seed"
          touch $tempfile
          for sim in $(seq $simspercfg); do
              echo "$configfile $seedcounter" >> $tempfile
              seedcounter=$((seedcounter+1))
          done
          mv $tempfile $seedfile
      done
  fi
  echo "Generating list of config paths and seeds... OK"
fi


# Take seeds from seedfiles, distribute them to config files, randomize and then split into .seed files
if [ -n "$configfiles" ] && [ -n  "$seedfiles" ] && [ -n "$shuffle" ]; then
    # Generate a master list with 2 columns: configs and seeds
    echo "Generating list config paths and given seeds..."
    superfile=$seedpath/superfile.txt
    touch $superfile
    for seedfile in $seedfiles; do
        seedbase=$(basename $seedfile .seed)
        configfile=$(find -L $configpath -type f -name $seedbase.cfg |  sort -g )
        if [ -z "$configfile" ]; then echo "Could not find a matching config file for seed file [ $seedfile ]. Searched with: find -L $configpath -type f -name $seedbase.cfg | Found: $match" ; exit 1; fi
        for seed in $(cat $seedfile); do
            echo "$configfile $seed" >> $superfile
        done
    done
    echo "Shuffling seeds ..."
    shuffledfile=$seedpath/shuffledfile.txt
    truncate -s 0 $shuffledfile
    cat $superfile | shuf --output $shuffledfile
    echo "Shuffling seeds ... OK"
    echo "Splitting list ..."
    seedcounter=$((startseed + 0))
    for configfile in $configfiles; do
        endseed=$seedcounter+$simspercfg
        seedfile="$seedpath/part.[$seedcounter-$endseed].seed"
        touch $seedfile
        # Move simspercfg lines from shuffledfile --> seedfile
        head -$simspercfg $shuffledfile > $seedfile && sed -i '1,$simspercfgd' $shuffledfile
        seedcounter=$((seedcounter+simspercfg))
    done
    rm $superfile $shuffledfile
    echo "Splitting list ... OK"
    echo "Generating list with config paths and given seeds... OK"
fi



# Now we have many .seed files that may need to be chunked into .job files of length simsperbatch
echo "Splitting seed files into job files..."
seedfiles=$(find -L $seedpath -type f -name '*.seed' |  sort -g )
for seedfile in $seedfiles; do
  seedbase=$(basename $seedfile .seed)
  if [ "$simspercfg" -gt "$simsperarray" ]; then
    split --lines=$simsperarray --additional-suffix=.job -d --suffix-length=3 $seedfile $seedpath/$seedbase.chunk-
    rm $seedfile
  else
    mv $seedfile $seedpath/$seedbase.job
  fi
done
echo "Splitting seed files into job files... OK"

# From this point on we are guaranteed to have a set of files in $seedpath/<file>.###.job
# Each .job file contains 2 columns with the corresponding config files and seeds
# The length of each .job file can vary!

jobfiles=$(find -L $seedpath -type f -name '*.job' |  sort -g )
if [ -z "$jobfiles" ] ; then
    echo "No simulation files found!"
    exit 1
fi


if [ -z "$dryrun" ]; then
  git log -1 >> job_report.txt
fi

for jobfile in $jobfiles; do

  numseeds=$(cat $jobfile | wc -l)

  if [ -n "$dryrun" ]; then
cat << EOF >&2
sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
--array=1-$numseeds:$simspertask \
run_jobarray.sh -e $exec -f $jobfile
EOF
  bash run_jobarray.sh -e $exec -f $jobfile -d
  else
    echo "sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
      --array=1-$numseeds:$simspertask \
      run_jobarray.sh -e $exec -f $jobfile" >> job_report.txt

    sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
      --array=1-$numseeds:$simspertask \
      run_jobarray.sh -e $exec -f $jobfile
  fi
done

