#!/bin/bash
PROGNAME=$0


usage() {
  cat << EOF >&2
Usage            : $PROGNAME [-option <argument>]

-a               : Choose microarchitecture for cxx and openblas. | core2 | nehalem | sandybridge | haswell | native | (default = haswell)
-b <build type>  : Build type: [ Release | RelWithDebInfo | Debug | Profile ]  (default = Release)
                   (Use the same build type you used with build.sh)
-c               : Clear CMake files before build (delete ./build)
-d               : Dry run
-g <compiler>    : Compiler        | GNU | Clang | (default = "")
-h               : Help. Shows this text.
-i <ON|OFF>      : Intel MKL use   | ON | OFF | (default = OFF)
-j <num_threads> : Number of threads used by CMake
-l <lib name>    : Clear library before build (i.e. delete ./libs/<lib name> and ./cmake-build-libs/<lib name>)
-L               : Clear all downloaded libraries before build (i.e. delete ./libs and ./cmake-build-libs)
-o <ON|OFF>      : OpenMP use      | ON | OFF | (default = OFF)
-s <ON|OFF>      : Shared libs     | ON | OFF | (default = OFF)
-t <target>      : DMRG++          | all | hdf5_test_target | arpack++_simple_test_target | arpack++_mps_test_target | (default = all)
-p <path>        : Path to gcc installation (default = )

EXAMPLE:
./build.sh -s OFF -a native -j 20 -o OFF -i ON -b Release -c -g Clang
EOF
  exit 1
}


target="all"
build="Release"
clear_cmake=""
clear_libs=""
march="haswell"
omp="OFF"
mkl="OFF"
shared="OFF"
compiler=""

while getopts a:b:cdg:hi:j:l:Lo:p:s:t: o; do
    case $o in
	    (a) march=$OPTARG;;
        (b) build=$OPTARG;;
        (c) clear_cmake="true";;
        (d) dryrun="true";;
        (g) compiler=$OPTARG;;
        (h) usage ;;
        (j) make_threads=$OPTARG;;
        (l) clear_lib+=("$OPTARG");;
        (L) clear_libs="true";;
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
    build_lowercase=$(echo $build | tr '[:upper:]' '[:lower:]')
	rm -rf ./libs-$build_lowercase ./build/$build/external-deps
else
    for lib in "${clear_lib[@]}"; do
        rm -r ./libs-$build_lowercase/$lib ./build/$build/external-deps/$lib
    done
fi






if [[ "$OSTYPE" == "linux-gnu" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "OS: Linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "OS: Mac OSX"
    echo "Checking if gcc-7 compiler is available"
    if brew ls gcc@7 | grep -q 'g++-7'; then
        echo " gcc-7 was found!"
        gccver=7
    elif brew ls gcc@8 | grep -q 'g++-8'; then
        echo " gcc-8 was found!"
        gccver=8
    else
        echo "Please install gcc (version 7 or 8) through brew."
        echo "Command:   brew install gcc@7"
    fi
    export CC=gcc-$gccver
    export CXX=g++-$gccver
    export FC=gfortran-$gccver
fi



if [[ "$HOSTNAME" == *"tetralith"* ]];then
    echo "Running on tetralith"
    conda activate dmrg
    module load buildenv-gcc/2018a-eb
    module load zlib
    module load CMake/3.15.2
#    module load GCCcore
    if [ "$compiler" = "GNU" ] ; then
        export CC=gcc
        export CXX=g++
    elif [ "$compiler" = "Clang" ] ; then
        module load Clang/8.0.0-GCCcore-8.2.0
        if [ -z "$gcc_toolchain" ] ; then
            gcc_toolchain=--gcc-toolchain=$EBROOTGCCCORE
        fi
        export CC=clang
        export CXX=clang++
    fi

elif [[ "$HOSTNAME" == *"raken"* ]];then
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


if [ -n "$dryrun" ]; then
    echo "Dry run build sequence"
else
    echo "Running build sequence"
fi

cat << EOF >&2
    cmake -E make_directory build/$build
    cd build/$build
    cmake -DCMAKE_BUILD_TYPE=$build -DMARCH=$march  -DUSE_OpenMP=$omp -DUSE_MKL=$mkl -DBUILD_SHARED_LIBS=$shared -DGCC_TOOLCHAIN=$gcc_toolchain  -G "CodeBlocks - Unix Makefiles" ../../
    cmake --build . --target $target -- -j $make_threads
EOF

if [ -z "$dryrun" ] ;then
    cmake -E make_directory build/$build
    cd build/$build
    cmake -DCMAKE_BUILD_TYPE=$build -DMARCH=$march  -DUSE_OpenMP=$omp -DUSE_MKL=$mkl -DBUILD_SHARED_LIBS=$shared -DGCC_TOOLCHAIN=$gcc_toolchain  -G "CodeBlocks - Unix Makefiles" ../../
    cmake --build . --target $target -- -j $make_threads
fi