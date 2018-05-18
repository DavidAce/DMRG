#!/bin/bash
PROGNAME=$0
default_mode="Release"
default_threads="2"
default_file=""
usage() {
  cat << EOF >&2
Usage: $PROGNAME [-t <target>][-m <mode>] [-j <num_threads>] [-i <input_file>]

-t <target>      : DMRG++    | all   | any test target
-m <mode>        : Release   | Debug
-j <num_threads> : Number of threads used by CMake
-i <input_file>  : Full or relative path to the input file
EOF
  exit 1
}

mode=$default_mode
threads=$default_threads
file=$default_file
while getopts t:m:j:i: o; do
  case $o in
    (t) target=$OPTARG;;
    (m) mode=$OPTARG;;
    (j) threads=$OPTARG;;
    (i) file=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"



./build/${mode}/DMRG++ $file
