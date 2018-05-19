#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-h ] [-t <target>] [-m <mode>] [-c]  [-l] [-j <num_threads>]
-c               : Clear CMake files before build (delete ./build)
-h               : Help. Shows this text.
-j <num_threads> : Number of threads used by CMake
-l               : Clear downloaded libraries (delete ./libs)
-m <mode>        : Release   | Debug
-t <target>      : DMRG++    | all   | any test target
EOF
  exit 1
}


target="all"
mode="Release"
clear_cmake=""
clear_libs=""
threads="2"

while getopts chj:lm:t: o; do
    case $o in
        (c) clear_cmake="true";;
        (h) usage ;;
        (j) threads=$OPTARG;;
        (l) clear_libs="true";;
        (m) mode=$OPTARG;;
        (t) target=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done
shift "$((OPTIND - 1))"


if [ "$clear_cmake" = "true" ]
then
    echo "Clearing CMake files from build."
	rm -rf ./build
fi

if [ "$clear_libs" = "true" ]
then
    echo "Clearing downloaded libraries."
	rm -rf ./libs
fi



if [[ "$OSTYPE" == "linux-gnu" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "OS: Linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "OS: Mac OSX"
    echo "Checking if gcc-7 compiler is available"
    if brew ls gcc@7 | grep -q 'g++-7'; then
        echo " gcc-7 was found!"
        if [ CC != "gcc-7" ]; then
            echo "Please export before running: "
            echo "  export CC=gcc-7"
            echo "  export CXX=g++-7"
            echo "  export FC=gfortran-7"
        fi
    else
        echo "Please install gcc (version 7 or higher) through brew."
        echo "Command:   brew install gcc@7"
#        exit 1
    fi
fi


echo "Starting Build"
echo "Target          :   $target"
echo "Build threads   :   $threads"
echo "Mode            :   $mode"

cmake -E make_directory build/$mode
cd build/$mode
cmake -DCMAKE_BUILD_TYPE=$mode -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target $target -- -j $threads
