#!/bin/bash
buildtype="Release"
dcmake_c_compiler=""
dcmake_cxx_compiler=""
if [[ "$@" == *"ebug"* ]]
then
	buildtype="Debug"
fi

if [[ "$@" == *"lean"* ]]
then
    echo "Cleaning build"
	rm -rf build
	rm -rf cmake/download_scripts/tmp/
    exit 0
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
    echo "Checking that GCC is installed in homebrew: [brew ls gcc | grep -q 'g++-7']"
    if brew ls gcc | grep -q 'g++-7'; then
        echo "GCC-7 is installed"
        export CC=gcc-7
        export CXX=g++-7
        dcmake_c_compiler="-DCMAKE_C_COMPILER=gcc-7"
        dcmake_cxx_compiler="-DCMAKE_CXX_COMPILER=g++-7"
    else
        echo "Please install GCC (version 7) through homebrew"
    fi
else
        echo "Could not identify OS"
fi


echo "Starting Build"
cmake -E make_directory build/${buildtype}
cd build/${buildtype}
cmake -DCMAKE_BUILD_TYPE=${buildtype} -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target DMRG++ -- -j 4
