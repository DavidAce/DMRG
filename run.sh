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

while getopts f:hm:t: o; do
      case $o in
        (f) file=$OPTARG;;
        (h) usage ;;
        (m) mode=$OPTARG;;
        (t) target=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage
  esac
done
shift "$((OPTIND - 1))"

echo "Running command:  ./build/$mode/$target $file"

./build/$mode/$target $file
