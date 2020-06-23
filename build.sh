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
-g | --compiler [=arg]          : Compiler        | GNU | Clang | Tau (default = "")
-G | --generator [=arg]         : CMake generator  | many options... | (default = "CodeBlocks - Unix Makefiles")
   | --gcc-toolchain [=arg]     : Path to GCC toolchain. Use with Clang if it can't find stdlib (defailt = none)
-h | --help                     : Help. Shows this text.
-j | --make-threads [=num]      : Number of threads used by Make build (default = 8)
-l | --clear-libs [=args]       : Clear libraries in comma separated list 'lib1,lib2...'. "all" deletes all.
-s | --enable-shared            : Enable shared library linking (default is static)
   | --enable-openmp            : Enable OpenMP
   | --enable-mkl               : Enable Intel MKL
   | --enable-lto               : Enable Link Time Optimization
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
                enable-lto\
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
enable_lto="OFF"
make_threads=8
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
    -a|--arch)                      march=$2                        ; echo " * Architecture            : $2"      ; shift 2 ;;
    -b|--build-type)                build_type=$2                   ; echo " * Build type              : $2"      ; shift 2 ;;
    -c|--clear-cmake)               clear_cmake="ON"                ; echo " * Clear CMake             : ON"      ; shift   ;;
    -l|--clear-libs)
            clear_libs=($(echo "$2" | tr ',' ' '))                  ; echo " * Clear libraries         : $2"      ; shift 2 ;;
    -d|--dry-run)                   dry_run="ON"                    ; echo " * Dry run                 : ON"      ; shift   ;;
       --download-method)           download_method=$2              ; echo " * Download method         : $2"      ; shift 2 ;;
    -f|--extra-flags)               extra_flags=$2                  ; echo " * Extra CMake flags       : $2"      ; shift 2 ;;
    -g|--compiler)                  compiler=$2                     ; echo " * C++ Compiler            : $2"      ; shift 2 ;;
    -G|--generator)                 generator=$2                    ; echo " * CMake generator         : $2"      ; shift 2 ;;
    -j|--make-threads)              make_threads=$2                 ; echo " * MAKE threads            : $2"      ; shift 2 ;;
    -s|--enable-shared)             enable_shared="ON"              ; echo " * Link shared libraries   : ON"      ; shift   ;;
       --enable-tests)              enable_tests="ON"               ; echo " * CTest Testing           : ON"      ; shift   ;;
    -t|--target)                    target=$2                       ; echo " * CMake Build target      : $2"      ; shift 2 ;;
       --enable-openmp)             enable_openmp="ON"              ; echo " * OpenMP                  : ON"      ; shift   ;;
       --enable-mkl)                enable_mkl="ON"                 ; echo " * Intel MKL               : ON"      ; shift   ;;
       --enable-lto)                enable_lto="ON"                 ; echo " * Link Time Optimization  : ON"      ; shift   ;;
       --no-modules)                no_modules="ON"                 ; echo " * Disable module load     : ON"      ; shift   ;;
       --prefer-conda)              prefer_conda="ON"               ; echo " * Prefer anaconda libs    : ON"      ; shift   ;;
    -v|--verbose)                   verbose="ON"                    ; echo " * Verbose makefiles       : ON"      ; shift   ;;
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


if [[ "$prefer_conda" =~ OFF|off|False|false ]] && [ -n "$CONDA_PREFIX" ] ; then
    if [ -f "$CONDA_PREFIX_1/etc/profile.d/conda.sh" ]; then
        source $CONDA_PREFIX_1/etc/profile.d/conda.sh
    fi
    if [ -f "$CONDA_PREFIX/etc/profile.d/conda.sh" ]; then
        source $CONDA_PREFIX/etc/profile.d/conda.sh
    fi
    counter=0
    while [ -n "$CONDA_PREFIX" ] && [  $counter -lt 10 ]; do
        let counter=counter+1
        echo "Deactivating conda environment $CONDA_PREFIX"
        conda deactivate
    done
    conda info --envs
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
        module load CMake/3.16.5
        module load buildenv-gcc/2018a-eb
        module load foss/2019a
        if [ "$enable_mkl" = "ON" ] ; then
            export MKLROOT=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a/mkl
            export EBROOTIMKL=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a
        elif [[ "$download_method" =~ find ]]; then
            module load OpenBLAS
        fi
        if [[ "$compiler" =~ Clang|clang|cl ]] ; then
            module try-load clang/6.0.1
        fi
        module list
    fi


    cmake --version
elif [[ "$HOSTNAME" == *"raken"* ]];then
    if [ -z "$no_module" ]; then
        module load CMake
        if [ "$enable_mkl" = "ON" ] ; then
            module load imkl
        fi
        if [[ "$download_method" =~ find ]] ; then
                module load HDF5/1.10.5-GCCcore-8.3.0
                if [ "$enable_mkl" = "OFF" ] ; then
                    module load OpenBLAS
                fi
        fi
        if [[ "$compiler" =~ Clang|clang|cl ]] ; then
            module load Clang
        fi
        module list
    fi
fi


if [[ "$compiler" =~ GCC|Gcc|gcc|cc|GNU|gnu|Gnu ]] ; then
    echo "Exporting compiler flags for GCC"
    export CC=gcc
    export CXX=g++
elif [[ "$compiler" =~ Clang|clang|cl ]] ; then
    echo "Exporting compiler flags for Clang"
    export CC=clang
    export CXX=clang++
elif [[ "$compiler" =~ Tau|tau ]] ; then
    echo "Exporting compiler flags for Tau"
    echo "Enabling shared linking for Tau"
    echo "Hint: If this doesnt work, try compiling normally with the plain compiler once first"
    enable_shared="ON"
    export PATH=/home/david/GitProjects/DMRG++/.tau/bin/ThinkStation/:$PATH
    export CC=tau_gcc
    export CXX=tau_g++
    export FC=tau_gfortran
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
          -DDMRG_ENABLE_MKL=$enable_mkl \
          -DDMRG_ENABLE_LTO=$enable_lto \
          $extra_flags
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
          -DDMRG_ENABLE_LTO=$enable_lto \
          $extra_flags \
           -G "$generator" \
           ../../
    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo ""
            echo "Exit code: $exit_code"
            echo "CMakeFiles/CMakeError.log:"
            echo ""
 #           cat CMakeFiles/CMakeError.log
            exit "$exit_code"
    fi

    if [ "$enable_tests" = "ON" ] ;then
        if [[ "$target" == *"test-"* ]]; then
            ctest --build-and-test ../../ tests --build-generator $generator --build-config $build_type --verbose  --build-target $target -R $target
        else
            ctest --build-and-test ../../ tests --build-generator $generator --build-config $build_type --verbose  --output-on-failure
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
#            cat CMakeFiles/CMakeError.log
            exit "$exit_code"
    fi
fi


