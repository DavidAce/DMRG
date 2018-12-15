#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-c] [-h ] [-j <num_threads>] [-l] [-m <mode>] [-t <target>] [-a <march>]

-a               : Choose microarchitecture for cxx and openblas. | core2 | nehalem | sandybridge | haswell | (default = haswell)
-c               : Clear CMake files before build (delete ./build)
-h               : Help. Shows this text.
-j <num_threads> : Number of threads used by CMake
-l               : Clear downloaded libraries before build (i.e. delete ./libs)
-m <mode>        : Release   | Debug | Profile |  (default = Release)
-t <target>      : DMRG++    | all   | any test target | (default = all)
EOF
  exit 1
}


target="all"
mode="Release"
clear_cmake=""
clear_libs=""
threads="8"
march="sandybridge"

while getopts a:chj:lm:t: o; do
    case $o in
	(a) march=$OPTARG;;
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
if [ $OPTIND -eq 1 ]; then echo "No flags were passed"; usage ;exit 1; fi
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
        if [[ "${CC}" != *"gcc-7"* ]]; then
            echo "Current CC is: ${CC}"
            echo "Please export before running: "
            echo "  export CC=gcc-7"
            echo "  export CXX=g++-7"
            echo "  export FC=gfortran-7"
        fi
    elif brew ls gcc@8 | grep -q 'g++-8'; then
        echo " gcc-8 was found!"
        if [[ "${CC}" != *"gcc-8"* ]]; then
            echo "Current CC is: ${CC}"
            echo "Please export before running: "
            echo "  export CC=gcc-8"
            echo "  export CXX=g++-8"
            echo "  export FC=gfortran-8"
        fi
    else
        echo "Please install gcc (version 7 or higher) through brew."
        echo "Command:   brew install gcc@7"
    fi
fi


echo "Starting Build"
echo "Micro arch.     :   $march"
echo "Target          :   $target"
echo "Build threads   :   $threads"
echo "Mode            :   $mode"

#module load CLANG_6.x.x
module load GNU_8.x.x
module load openblas_${march}_v0.3.3
module load arpack-ng_${march}_3.6.2
module load armadillo-9.200.x
module load arpack++
module load hdf5_1.10.3
module load gsl_2.4
module load eigen3_3.3.5

cmake -E make_directory build/$mode
cd build/$mode
cmake -DCMAKE_BUILD_TYPE=$mode -DMARCH=$march  -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target $target -- -j $threads
