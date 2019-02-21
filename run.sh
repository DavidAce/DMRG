#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-f <input_file>] [-h] [-m <mode>] [-t <target>]

-f <input_file>  : Full or relative path to the input file
-h               : Help. Shows this text.
-m <mode>        : Release   | Debug | (default = Release) -- (Use the same mode you used with build.sh)
-t <target>      : DMRG++    | all   | any test target | (default = DMRG++)
EOF
  exit 1
}

target="DMRG++"
mode="Release"
file="input.cfg"
threads="8"

while getopts f:hm:t: o; do
      case $o in
        (f) file=$OPTARG;;
        (h) usage ;;
        (m) mode=$OPTARG;;
        (t) target=$OPTARG;;
        (n) threads=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        (*) usage ;;
      esac
done

if [ $OPTIND -eq 1 ]; then echo "No flags were passed"; usage ;exit 1; fi





if [[ "$HOSTNAME" == *"tetralith"* ]];then
     sbatch run_triolith.sh ./build/$mode/$target $file
#    ./run_tetralith.sh ./build/$mode/$target $file
else
    echo "Threads         :   $threads"
    export OMP_NUM_THREADS=$threads

    shift "$((OPTIND - 1))"
    ulimit -c unlimited
    echo "Running command:  ./build/$mode/$target $file"

    ./build/$mode/$target $file
fi

