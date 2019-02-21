#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-c] [-h ] [-j <num_threads>] [-l] [-m <mode>] [-t <target>] [-a <march>]

-a               : Choose microarchitecture for cxx and openblas. | core2 | nehalem | sandybridge | haswell | native | (default = sandybridge)
-c               : Clear CMake files before build (delete ./build)
-g <compiler>    : Compiler        | GNU | Clang | (default = "")
-h               : Help. Shows this text.
-i <ON|OFF>      : Intel MKL use   | ON | OFF | (default = OFF)
-j <num_threads> : Number of threads used by CMake
-l               : Clear downloaded libraries before build (i.e. delete ./libs)
-m <mode>        : Release         | Debug | Profile |  (default = Release)
-o <ON|OFF>      : OpenMP use      | ON | OFF | (default = OFF)
-s <ON|OFF>      : Static linking  | ON | OFF | (default = ON)
-t <target>      : DMRG++          | all | hdf5_test_target | arpack++_simple_test_target | arpack++_mps_test_target | (default = all)
EOF
  exit 1
}


target="all"
mode="Release"
clear_cmake=""
clear_libs=""
make_threads="8"
march="sandybridge"
omp="OFF"
mkl="OFF"
static="ON"
compiler=""

while getopts a:cg:hi:j:lm:o:s:t: o; do
    case $o in
	    (a) march=$OPTARG;;
        (c) clear_cmake="true";;
        (g) compiler=$OPTARG;;
        (h) usage ;;
        (j) make_threads=$OPTARG;;
        (l) clear_libs="true";;
        (m) mode=$OPTARG;;
        (t) target=$OPTARG;;
        (o) omp=$OPTARG;;
        (i) mkl=$OPTARG;;
        (s) static=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done

if [ $OPTIND -eq 1 ] ; then
    echo "No flags were passed"; usage ;exit 1;
fi
shift "$((OPTIND - 1))"


if [ "$clear_cmake" = true ] ; then
    echo "Clearing CMake files from build."
	rm -rf ./build
fi

if [ "$clear_libs" = true ] ; then
    echo "Clearing downloaded libraries."
	rm -rf ./libs
fi


if [ "$compiler" = "GNU" ] ; then
    module load GNU_8.x.x
elif [ "$compiler" = "Clang" ] ; then
    module load CLANG_7
fi


if [ "$mkl" = "ON" ] ; then
    module load intel-mkl-2019.1
else
    module load openblas_${march}_v0.3.4
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



if [[ "$hostname" == *"tetralith"* ]];then
    echo "Running on tetralith"
    module load CMake/3.12.1
    module load buildenv-gcc/7.3.0-bare
else
    module load arpack-ng_${march}_3.6.2
    module load arpack++
    module load hdf5_1.10.3
    module load gsl_2.4
    module load eigen3_3.3.5
fi




echo "Starting Build"
echo "Compiler        :   $compiler"
echo "Micro arch.     :   $march"
echo "Target          :   $target"
echo "Build threads   :   $make_threads"
echo "Mode            :   $mode"
echo "OpenMP          :   $omp"
echo "Intel MKL       :   $mkl"
echo "Static build    :   $static"



#module load CLANG_6.x.x


cmake -E make_directory build/$mode
cd build/$mode
cmake -DCMAKE_BUILD_TYPE=$mode -DMARCH=$march  -DUSE_OpenMP=$omp -DUSE_MKL=$mkl -DSTATIC_BUILD=$static  -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target $target -- -j $make_threads
