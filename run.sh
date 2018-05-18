#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-t <target>][-m <mode>] [-f <input_file>]

-t <target>      : DMRG++    | all   | any test target
-m <mode>        : Release   | Debug
-f <input_file>  : Full or relative path to the input file
EOF
  exit 1
}

target="all"
mode="Release"
file="input.cfg"

while getopts t:m:f: o; do
  case $o in
    (t) target=$OPTARG;;
    (m) mode=$OPTARG;;
    (f) file=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"



./build/$mode/$target $file
