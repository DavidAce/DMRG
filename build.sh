#!/bin/bash
buildtype="Release"
dcmake_c_compiler=""
dcmake_cxx_compiler=""
dcmake_fortran_compiler=""

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
    echo "Checking if fortran compiler is available"
    fortran_compiler=$(brew ls gcc | grep -e '/bin/gfortran-[0-9]' | head -n 1)
    if [ -z "${fortran_compiler}" ];then
        echo "Fortran compiler is not available. Automatic library installation (e.g. GSL and Arpack) will not work"
        echo "Either install gfortran (bundled with gcc) or install these libraries manually."
    fi

    echo "Checking if libstdc++ is available"
    standard_library=$(brew ls gcc | grep -e '[0-9]/lib/gcc/[0-9]/libstdc++.a' | head -n 1)
    if [ -z "${standard_library}" ];then
        echo "Please install standard library available in gcc first."
        echo "  brew install gcc"
        exit
    fi

    echo "Checking if clang is installed through brew..."
    if brew ls llvm | grep -q 'clang++'; then
        c_compiler=$(brew ls llvm | grep -e '/bin/clang' | head -n 1)
        cxx_compiler=$(brew ls llvm | grep -e '/bin/clang++' | head -n 1)
        echo "clang is installed"

    elif brew ls gcc | grep -q 'g++'; then
        c_compiler=$(brew ls gcc | grep -e '/bin/gcc-[0-9]' | head -n 1)
        cxx_compiler=$(brew ls gcc | grep  -e '/bin/g++-[0-9]' | head -n 1)
        echo "gcc is installed"
    else
        echo "Please install gcc (version 7 or higher) or gcc + llvm (version 6 or higher) through brew"
        echo "  brew install gcc"
        echo "  brew install llvm"
        exit
    fi
    export LIBRARY_PATH=${standard_library}
    export CC=${c_compiler}
    export CXX=${cxx_compiler}
    export FC=${fortran_compiler}
    dcmake_c_compiler="-DCMAKE_C_COMPILER=${c_compiler}"
    dcmake_cxx_compiler="-DCMAKE_CXX_COMPILER=${cxx_compiler}"
    dcmake_fortran_compiler="-DCMAKE_Fortran_COMPILER=${fortran_compiler}"
    dcmake_library_path="-DCMAKE_LIBRARY_PATH=${standard_library}"
else
        echo "Could not identify OS"
        exit
fi


echo "Starting Build"
cmake -E make_directory build/${buildtype}
cd build/${buildtype}
cmake  ${dcmake_c_compiler} ${dcmake_cxx_compiler} ${dcmake_fortran_compiler} ${dcmake_library_path} -DCMAKE_BUILD_TYPE=${buildtype} -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target DMRG++ -- -j 4
