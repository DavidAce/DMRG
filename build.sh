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
-l               : Clear downloaded libraries before build (i.e. delete ./libs and ./cmake-build-libs)
-m <mode>        : Release         | Debug | Profile |  (default = Release)
-o <ON|OFF>      : OpenMP use      | ON | OFF | (default = OFF)
-s <ON|OFF>      : Shared libs     | ON | OFF | (default = OFF)
-t <target>      : DMRG++          | all | hdf5_test_target | arpack++_simple_test_target | arpack++_mps_test_target | (default = all)
-w <path>        : Path to gcc installation (default = "")
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
shared="OFF"
compiler=""
gcc_toolchain=""

while getopts a:cg:hi:j:lm:o:s:t: o; do
    case $o in
	    (a) march=$OPTARG;;
        (c) clear_cmake="true";;
        (g) compiler=$OPTARG;;
        (h) usage ;;
        (j) make_threads=$OPTARG;;
        (l) clear_libs="true";;
        (m) mode=$OPTARG;;
        (o) omp=$OPTARG;;
        (i) mkl=$OPTARG;;
        (s) shared=$OPTARG;;
        (t) target=$OPTARG;;
        (w) gcc_toolchain=$OPTARG;;
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
	rm -rf ./libs ./cmake-build-libs
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



if [[ "$HOSTNAME" == *"tetralith"* ]];then
    echo "Running on tetralith"
    module load buildenv-gcc/2018a-eb
    module load zlib/1.2.8
    if [ -z "$gcc_toolchain" ] ; then
        gcc_toolchain=/software/sse/easybuild/prefix/software/GCCcore/7.3.0
    fi
    module load clang/6.0.1
    module load GCCcore/7.3.0
    module load CMake/3.12.1
    source activate dmrg
    #export CC=gcc
    #export CXX=g++
    export CC=clang
    export CXX=clang++



elif [[ "$HOSTNAME" == *"anderson"* ]];then
    module load CMake
    if [ "$mkl" = "ON" ] ; then
        module load imkl
    else
        module load OpenBLAS
    fi
    module load arpack-ng
    module load ARPACK++
    module load HDF5/1.10.5-GCCcore-8.2.0
    module load Eigen
    module load glog

    if [ "$compiler" = "GNU" ] ; then
        module load GCCcore
    elif [ "$compiler" = "Clang" ] ; then
        module load Clang
    fi
fi




echo "Starting Build"
echo "Compiler        :   $compiler"
echo "CC              :   $CC"
echo "CXX             :   $CXX"
echo "Micro arch.     :   $march"
echo "Target          :   $target"
echo "Build threads   :   $make_threads"
echo "Mode            :   $mode"
echo "OpenMP          :   $omp"
echo "Intel MKL       :   $mkl"
echo "Shared build    :   $shared"
echo "gcc toolchain   :   $gcc_toolchain"
echo "CMake version   :   $(cmake --version) at $(which cmake)"



cmake -E make_directory build/$mode
cd build/$mode
cmake -DCMAKE_BUILD_TYPE=$mode -DMARCH=$march  -DUSE_OpenMP=$omp -DUSE_MKL=$mkl -DBUILD_SHARED_LIBS=$shared -DGCC_TOOLCHAIN=$gcc_toolchain  -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target $target -- -j $make_threads
