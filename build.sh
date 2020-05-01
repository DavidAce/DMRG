#!/bin/bash

PROGNAME=$0

usage() {
  cat << EOF >&2
Usage            : $PROGNAME [-option | --option ] <=argument>

-a | --arch [=arg]              : Choose microarchitecture | core2 | nehalem | sandybridge | haswell | native | (default = haswell)
-b | --build-type [=arg]        : Build type: [ Release | RelWithDebInfo | Debug | Profile ]  (default = Release)
-c | --clear-cmake              : Clear CMake files before build (delete ./build)
-d | --dry-run                  : Dry run
   | --download-method          : Download method for dependencies [ find | fetch | find-or-fetch | conan ] (default = find)
-f | --extra-flags [=arg]       : Extra CMake flags (defailt = none)
-g | --compiler [=arg]          : Compiler        | GNU | Clang | (default = "")
-G | --generator [=arg]         : CMake generator  | many options... | (default = "CodeBlocks - Unix Makefiles")
   | --gcc-toolchain [=arg]     : Path to GCC toolchain. Use with Clang if it can't find stdlib (defailt = none)
-h | --help                     : Help. Shows this text.
-j | --make-threads [=num]      : Number of threads used by Make build (default = 8)
-l | --clear-libs [=args]       : Clear libraries in comma separated list 'lib1,lib2...'. "all" deletes all.
-s | --enable-shared            : Enable shared library linking (default is static)
   | --enable-openmp            : Enable OpenMP
   | --enable-mkl               : Enable Intel MKL
-t | --target [=args]           : Select build target [ CMakeTemplate | all-tests | test-<name> ]  (default = none)
   | --enable-tests             : Enable CTest tests
   | --prefer-conda             : Prefer libraries from anaconda
   | --no-modules               : Disable use of "module load"
-v | --verbose                  : Verbose makefiles
EXAMPLE:
./build.sh --arch native -b Release  --make-threads 8   --enable-shared  --with-openmp --with-eigen3  --download-method=find
EOF
  exit 1
}


# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"   -o ha:b:cl:df:g:G:j:st:v \
                --long "\
                help\
                arch:\
                build-type:\
                target:\
                clear-cmake\
                clear-libs:\
                compiler:\
                dry-run\
                download-method:\
                enable-tests\
                enable-shared\
                gcc-toolchain:\
                make-threads:\
                enable-openmp\
                enable-mkl\
                no-modules\
                prefer-conda\
                verbose\
                generator\
                extra-flags:\
                "  -- "$@")

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ]; then exit 1 ; fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

build_type="Release"
target="all"
march="haswell"
enable_shared="OFF"
download_method="find"
enable_tests="OFF"
enable_openmp="OFF"
enable_mkl="OFF"
make_threads=1
prefer_conda="OFF"
verbose="OFF"
generator="CodeBlocks - Unix Makefiles"
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
#$1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
echo "Enabled options:"
while true;
do
  case "$1" in
    -h|--help)                      usage                                                                          ; shift   ;;
    -a|--arch)                      march=$2                        ; echo " * Architecture             : $2"      ; shift 2 ;;
    -b|--build-type)                build_type=$2                   ; echo " * Build type               : $2"      ; shift 2 ;;
    -c|--clear-cmake)               clear_cmake="ON"                ; echo " * Clear CMake              : ON"      ; shift   ;;
    -l|--clear-libs)
            clear_libs=($(echo "$2" | tr ',' ' '))                  ; echo " * Clear libraries          : $2"      ; shift 2 ;;
    -d|--dry-run)                   dry_run="ON"                    ; echo " * Dry run                  : ON"      ; shift   ;;
       --download-method)           download_method=$2              ; echo " * Download method          : $2"      ; shift 2 ;;
    -f|--extra-flags)               extra_flags=$2                  ; echo " * Extra CMake flags        : $2"      ; shift 2 ;;
    -g|--compiler)                  compiler=$2                     ; echo " * C++ Compiler             : $2"      ; shift 2 ;;
    -G|--generator)                 generator=$2                    ; echo " * CMake generator          : $2"      ; shift 2 ;;
       --gcc-toolchain)             gcc_toolchain=$2                ; echo " * GCC toolchain            : $2"      ; shift 2 ;;
    -j|--make-threads)              make_threads=$2                 ; echo " * MAKE threads             : $2"      ; shift 2 ;;
    -s|--enable-shared)             enable_shared="ON"              ; echo " * Link shared libraries    : ON"      ; shift   ;;
       --enable-tests)              enable_tests="ON"               ; echo " * CTest Testing            : ON"      ; shift   ;;
    -t|--target)                    target=$2                       ; echo " * CMake Build target       : $2"      ; shift 2 ;;
       --enable-openmp)             enable_openmp="ON"              ; echo " * Enable OpenMP            : ON"      ; shift   ;;
       --enable-mkl)                enable_mkl="ON"                 ; echo " * Enable Intel enable_mkl  : ON"      ; shift   ;;
       --no-modules)                no_modules="ON"                 ; echo " * Disable module load      : ON"      ; shift   ;;
       --prefer-conda)              prefer_conda="ON"               ; echo " * Prefer anaconda libs:    : ON"      ; shift   ;;
    -v|--verbose)                   verbose="ON"                    ; echo " * Verbose makefiles        : ON"      ; shift   ;;
    --) shift; break;;
  esac
done


if [ $OPTIND -eq 0 ] ; then
    echo "No flags were passed"; usage ;exit 1;
fi
shift "$((OPTIND - 1))"


if  [ -n "$clear_cmake" ] ; then
    echo "Clearing CMake files from build."
	rm -rf ./build/$build_type/CMakeCache.txt
fi

build_type_lower=$(echo $build_type | tr '[:upper:]' '[:lower:]')
for lib in "${clear_libs[@]}"; do
    if [[ "$lib" == "all" ]]; then
        echo "Clearing all installed libraries"
        rm -r ./build/$build_type/dmrg-deps-build/*
        rm -r ./build/$build_type/dmrg-deps-install/*
    else
        echo "Clearing library: $lib"
        rm -r ./build/$build_type/dmrg-deps-build/$lib
        rm -r ./build/$build_type/dmrg-deps-install/$lib
    fi
done

if [[ ! "$download_method" =~ find|fetch|conan ]]; then
    echo "Download method unsupported: $download_method"
    exit 1
fi




if [[ "$OSTYPE" == "darwin"* ]]; then
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
    if [ -z "$no_module" ]; then
        module load zlib
        module load GCC/8.2.0-2.31.1
        if [ "$compiler" = "Clang" ] ; then
            module load Clang/8.0.0-GCCcore-8.2.0
            if [ -z "$gcc_toolchain" ] ; then gcc_toolchain=--gcc-toolchain=$EBROOTGCCCORE ; fi
        fi
    fi

    if [ "$compiler" = "GNU" ] ; then
        export CC=gcc
        export CXX=g++
    elif [ "$compiler" = "Clang" ] ; then
        export CC=clang
        export CXX=clang++
    fi

elif [[ "$HOSTNAME" == *"raken"* ]];then
    if [ -z "$no_module" ]; then
        if [ "$enable_mkl" = "ON" ] ; then module load imkl; else module load OpenBLAS; fi
        module load arpack-ng
        module load ARPACK++
        module load HDF5/1.10.5-GCCcore-8.2.0
        module load Eigen # We want our own patched eigen though.
        module load CMake
        module load GCCcore
        if [ "$compiler" = "Clang" ] ; then
            module load Clang
            if [ -z "$gcc_toolchain" ] ; then gcc_toolchain=--gcc-toolchain=$EBROOTGCCCORE ; fi
        fi
        module list
    fi

    if [ "$compiler" = "GNU" ] ; then
        export CC=gcc
        export CXX=g++
    elif [ "$compiler" = "Clang" ] ; then
        export CC=clang
        export CXX=clang++
    fi
fi

export MAKEFLAGS=-j$make_threads
export CMAKE_BUILD_PARALLEL_LEVEL=$make_threads



echo " * Compiler                 :   $compiler"
echo " * MAKE threads (ext builds):   $make_threads"
echo " * CC                       :   $CC"
echo " * CXX                      :   $CXX"
echo " * CMake version            :   $(cmake --version) at $(which cmake)"

if [ -n "$dry_run" ]; then
    echo "Dry run build sequence"
else
    echo "Running build sequence"
fi

cat << EOF >&2
    cmake -E make_directory build/$build_type
    cd build/$build_type
    cmake -DCMAKE_BUILD_TYPE=$build_type
          -DBUILD_SHARED_LIBS=$enable_shared
          -DCMAKE_VERBOSE_MAKEFILE=$verbose
          -DDMRG_PRINT_INFO=$verbose
          -DDMRG_DOWNLOAD_METHOD=$download_method
          -DDMRG_PREFER_CONDA_LIBS:BOOL=$prefer_conda
          -DDMRG_MARCH=$march
          -DDMRG_ENABLE_TESTS:BOOL=$enable_tests
          -DDMRG_ENABLE_OPENMP=$enable_openmp
          -DDMRG_ENABLE_MKL=$enable_mkl
          -DGCC_TOOLCHAIN=$gcc_toolchain
           -G $generator
           ../../
    cmake --build . --target $target --parallel $make_threads
EOF

if [ -z "$dry_run" ] ;then
    cmake -E make_directory build/$build_type
    cd build/$build_type
    cmake -DCMAKE_BUILD_TYPE=$build_type \
          -DBUILD_SHARED_LIBS=$enable_shared \
          -DCMAKE_VERBOSE_MAKEFILE=$verbose \
          -DDMRG_PRINT_INFO=$verbose \
          -DDMRG_PRINT_CHECKS=$verbose \
          -DDMRG_DOWNLOAD_METHOD=$download_method \
          -DDMRG_PREFER_CONDA_LIBS:BOOL=$prefer_conda \
          -DDMRG_MARCH=$march \
          -DDMRG_ENABLE_TESTS:BOOL=$enable_tests \
          -DDMRG_ENABLE_OPENMP=$enable_openmp \
          -DDMRG_ENABLE_MKL=$enable_mkl \
          -DGCC_TOOLCHAIN=$gcc_toolchain \
           -G "$generator" \
           ../../
    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo ""
            echo "Exit code: $exit_code"
            echo "CMakeFiles/CMakeError.log:"
            echo ""
            cat CMakeFiles/CMakeError.log
            exit "$exit_code"
    fi

    if [ "$enable_tests" = "ON" ] ;then
        if [[ "$target" == *"test-"* ]]; then
            ctest  --build-config $build_type --verbose --build-and-test --build-target $target -R $target
        else
           ctest  --build-config $build_type --build-and-test --output-on-failure
        fi
    fi

    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo "Exit code: $exit_code"
            exit "$exit_code"
    fi


    cmake --build . --target $target --parallel $make_threads
    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo ""
            echo "Exit code: $exit_code"
            echo "CMakeFiles/CMakeError.log:"
            echo ""
            cat CMakeFiles/CMakeError.log
            exit "$exit_code"
    fi
fi


