#!/bin/bash 

buildtype="Release"
if [[ "$@" == *"ebug"* ]] 
then
	buildtype="Debug" 
fi
if [[ "$@" == *"lean"* ]] 
then
	echo "Cleaning build"
	rm -rf build
	exit 0 
fi 
if [[ "$HOSTNAME" == *"triolith"* ]] 
then
    echo "We're on triolith!";
    module add cmake/3.6.1
    module load_MPS buildenv-intel/2016-3
    export CC=/software/apps/gcc/5.3.0/build01/bin/gcc
    export CXX=/software/apps/gcc/5.3.0/build01/bin/g++ 
else
    echo "We're on my pc!" 
fi 
mkdir build 
cd build 
mkdir ${buildtype} 
cd ${buildtype} 
echo "Starting Build" 
cmake -DCMAKE_BUILD_TYPE=${buildtype} ../../
make
