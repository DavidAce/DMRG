#!/bin/bash
PROGNAME=$0
dcmake_c_compiler=""
dcmake_cxx_compiler=""
default_target="all"
default_mode="Release"
default_threads="2"

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-t <target>] [-m <mode>] [-c]  [-l] [-j <num_threads>]

-t <target>      : DMRG++    | all   | any test target
-m <mode>        : Release   | Debug
-c               : Clear CMake files before build (delete ./build)
-l               : Clear downloaded libraries (delete ./libs)
-j <num_threads> : Number of threads used by CMake
EOF
  exit 1
}

target=$default_target
mode=$default_mode
clear_cmake=""
clear_libs=""
threads=$default_threads
while getopts t:m:clj: o; do
  case $o in
    (t) target=$OPTARG;;
    (m) mode=$OPTARG;;
    (c) clear_cmake="true";;
    (l) clear_libs="true";;
    (j) threads=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"


if [ "${clear_cmake}" = "true" ]
then
    echo "Clearing CMake files from build."
	rm -rf ./build
fi

if [ "${clear_libs}" = "true" ]
then
    echo "Clearing downloaded libraries."
	rm -rf ./libs
fi

if [[ "$OSTYPE" == "linux-gnu" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "OS: Linux"
    if [[ "$HOSTNAME" == *"triolith"* ]]; then
        echo "We're on triolith!";
        module add cmake/3.6.1
        module load buildenv-intel/2016-3
        export CC=/software/apps/gcc/5.3.0/build01/bin/gcc
        export CXX=/software/apps/gcc/5.3.0/build01/bin/g++
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "OS: Mac OSX"
    echo "Checking if gcc compiler is available"

    if brew ls gcc | grep -q 'g++'; then
        c_compiler=$(brew ls gcc | grep -e '/bin/gcc-[0-9]' | head -n 1)
        cxx_compiler=$(brew ls gcc | grep  -e '/bin/g++-[0-9]' | head -n 1)
        fortran_compiler=$(brew ls gcc | grep -e '/bin/gfortran-[0-9]' | head -n 1)
        echo " -- gcc was found!"


#    echo "Checking if libstdc++ is available"
#    standard_library=$(brew ls gcc | grep -e '[0-9]/lib/gcc/[0-9]/libstdc++.a' | head -n 1)
#    standard_library_path=$(dirname ${standard_library})
#    if [ -z "${standard_library}" ];then
#        echo "Please install standard library available in gcc first."
#        echo "  brew install gcc"
#        exit
#    fi
#
#    if [ -z "${fortran_compiler}" ];then
#        echo "Fortran compiler is not available. Automatic library installation (e.g. GSL and Arpack) will not work"
#        echo "Either install gfortran (bundled with gcc) or install these libraries manually."
#    fi
#    echo "Checking if clang is installed through brew..."
#    if brew ls llvm | grep -q 'clang++'; then
#        c_compiler=$(brew ls llvm | grep -e '/bin/clang' | head -n 1)
#        cxx_compiler=$(brew ls llvm | grep -e '/bin/clang++' | head -n 1)
#        echo "clang is installed"


    else
        echo "Please install gcc (version 7 or higher) through brew"
        echo "  brew install gcc"
        exit 1
    fi
    export CC=${c_compiler}
    export CXX=${cxx_compiler}
    export FC=${fortran_compiler}
    dcmake_c_compiler="-DCMAKE_C_COMPILER=${c_compiler}"
    dcmake_cxx_compiler="-DCMAKE_CXX_COMPILER=${cxx_compiler}"
    dcmake_fortran_compiler="-DCMAKE_Fortran_COMPILER=${fortran_compiler}"
else
        echo "Could not identify OS"
        exit 1
fi


echo "Starting Build"
cmake -E make_directory build/$mode
cd build/$mode
cmake $dcmake_c_compiler $dcmake_cxx_compiler -DCMAKE_BUILD_TYPE=$mode -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target ${target} -- -j ${threads}
