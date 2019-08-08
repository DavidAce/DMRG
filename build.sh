#!/bin/bash
PROGNAME=$0


usage() {
  cat << EOF >&2
Usage            : $PROGNAME [-option <argument>]

-a               : Choose microarchitecture for cxx and openblas. | core2 | nehalem | sandybridge | haswell | native | (default = sandybridge)
-b <build type>  : Build type | Release | RelWithDebInfo | Debug | Profile |  (default = Release)
-c               : Clear CMake files before build (delete ./build)
-g <compiler>    : Compiler        | GNU | Clang | (default = "")
-h               : Help. Shows this text.
-i <ON|OFF>      : Intel MKL use   | ON | OFF | (default = OFF)
-j <num_threads> : Number of threads used by CMake
-l               : Clear downloaded libraries before build (i.e. delete ./libs and ./cmake-build-libs)
-o <ON|OFF>      : OpenMP use      | ON | OFF | (default = OFF)
-s <ON|OFF>      : Shared libs     | ON | OFF | (default = OFF)
-t <target>      : DMRG++          | all | hdf5_test_target | arpack++_simple_test_target | arpack++_mps_test_target | (default = all)
-p <path>        : Path to gcc installation (default = )

EXAMPLE:
./build.sh -s OFF -a native -j 20 -o OFF  -i ON -b Release -c  -g Clang
EOF
  exit 1
}


target="all"
build="Release"
clear_cmake=""
clear_libs=""
march="sandybridge"
omp="OFF"
mkl="OFF"
shared="OFF"
compiler=""

while getopts a:b:cg:hi:j:lo:p:s:t: o; do
    case $o in
	    (a) march=$OPTARG;;
        (b) build=$OPTARG;;
        (c) clear_cmake="true";;
        (g) compiler=$OPTARG;;
        (h) usage ;;
        (j) make_threads=$OPTARG;;
        (l) clear_libs="true";;
        (o) omp=$OPTARG;;
        (p) gcc_toolchain=--gcc-toolchain=$OPTARG;;
        (i) mkl=$OPTARG;;
        (s) shared=$OPTARG;;
        (t) target=$OPTARG;;
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
    conda activate dmrg
    module load buildenv-gcc/2018a-eb
    module load zlib
    module load CMake/3.12.1
    module load GCCcore
    if [ "$compiler" = "GNU" ] ; then
        export CC=gcc
        export CXX=g++
    elif [ "$compiler" = "Clang" ] ; then
        module load Clang
        if [ -z "$gcc_toolchain" ] ; then
            gcc_toolchain=--gcc-toolchain=$EBROOTGCCCORE
        fi
        export CC=clang
        export CXX=clang++
    fi



elif [[ "$HOSTNAME" == *"anderson"* ]];then
    module load CMake
    if [ "$mkl" = "ON" ] ; then
        module load imkl
    else
        module load OpenBLAS
    fi
    module load XZ/5.2.4-GCCcore-8.2.0 
    module load arpack-ng
    module load ARPACK++
    module load HDF5/1.10.5-GCCcore-8.2.0
    module load Eigen
    module load gflags
    module load glog
    module load CMake
    module load GCCcore
    module list
    if [ "$compiler" = "GNU" ] ; then
        export CC=gcc
        export CXX=g++
    elif [ "$compiler" = "Clang" ] ; then
        module load Clang
        if [ -z "$gcc_toolchain" ] ; then
            gcc_toolchain=--gcc-toolchain=$EBROOTGCCCORE
        fi
        export CC=clang
        export CXX=clang++
    fi
fi




echo "Starting Build"
echo "Compiler        :   $compiler"
echo "CC              :   $CC"
echo "CXX             :   $CXX"
echo "Micro arch.     :   $march"
echo "Target          :   $target"
echo "Build threads   :   $make_threads"
echo "Build Type      :   $build"
echo "OpenMP          :   $omp"
echo "Intel MKL       :   $mkl"
echo "Shared build    :   $shared"
echo "gcc toolchain   :   $gcc_toolchain"
echo "CMake version   :   $(cmake --version) at $(which cmake)"



cmake -E make_directory build/$build
cd build/$build
cmake -DCMAKE_BUILD_TYPE=$build -DMARCH=$march  -DUSE_OpenMP=$omp -DUSE_MKL=$mkl -DBUILD_SHARED_LIBS=$shared -DGCC_TOOLCHAIN=$gcc_toolchain  -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target $target -- -j $make_threads
