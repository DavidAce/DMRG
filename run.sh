#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-f <input_file>] [-h] [-m <mode>] [-t <target>]
-b <build type>  : Build type: [ Release | RelWithDebInfo | Debug | Profile ]  (default = Release)
                   (Use the same build type you used with build.sh)
-d               : Dry run
-e <execname>    : DMRG++    | all   | any test target | (default = DMRG++)
-f <input_file>  : Full or relative path to the input file
-h               : Help. Shows this text
-n <num threads> : Number of OpenMP threads
-r <seed>        : Seed for random number generator used for the model (default = 0 )
-s <seed>        : Seed for random number generator used for initial state (default = 0 )
EOF
  exit 1
}

execname="DMRG++"
build="Release"
file="input.cfg"
threads="8"
modelseed=0
stateseed=0
while getopts b:de:f:hn:r:s: o; do
      case $o in
        (b) build=$OPTARG;;
        (d) dryrun=true;;
        (e) execname=$OPTARG;;
        (f) file=$OPTARG;;
        (h) usage ;;
        (n) threads=$OPTARG;;
        (r) modelseed=$OPTARG;;
        (s) stateseed=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        (*) usage ;;
      esac
done

if [ $OPTIND -eq 1 ]; then echo "No flags were passed"; usage ;exit 1; fi


shift "$((OPTIND - 1))"
ulimit -c unlimited



echo "Threads : $threads"
export OMP_NUM_THREADS=$threads

exec=../build/$build/$execname
if [ -f "$exec" ]; then
    echo "Found executable: $exec"
else
    echo "Executable does not exist: $exec"
    exit 1
fi


if [ -n "$dryrun" ] ; then
    echo "Dry run command:  $exec -i $file -r $modelseed -s $stateseed"
else
    echo "Running command:  $exec -i $file -r $modelseed -s $stateseed"
    $exec -i $file -r $modelseed -s $stateseed
fi


