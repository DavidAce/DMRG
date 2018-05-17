#!/bin/bash
PROGNAME=$0
default_mode="Release"
default_threads="2"
default_file=""
usage() {
  cat << EOF >&2
Usage: $PROGNAME [-m <mode>] [-j <num_threads>] [-i <input_file>]

-m <mode>        : Release   | Debug
-j <num_threads> : Number of threads used by CMake
-i <input_file>  : Full or relative path to the input file
EOF
  exit 1
}

mode=$default_mode
threads=$default_threads
file=$default_file
while getopts m:j:i: o; do
  case $o in
    (m) mode=$OPTARG;;
    (j) threads=$OPTARG;;
    (i) file=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"



./build/${mode}/DMRG++ $file
