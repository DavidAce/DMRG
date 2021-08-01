#!/bin/bash

PROGNAME=$0

usage() {
  cat << EOF >&2
Usage            : $PROGNAME [-option | --option ] <=argument>

-a | --arch [=arg]              : Choose microarchitecture | core2 | nehalem | sandybridge | haswell | native | (default = haswell)
-b | --build-type [=arg]        : Build type: [ Release | RelWithDebInfo | Debug | Profile ]  (default = Release)
-c | --clear-cmake              : Clear CMake files before build (delete ./build)
-d | --dry-run                  : Dry run
   | --default-tetralith        : Set all settings for building on Tetralith
   | --default-kraken           : Set all settings for building on Kraken
   | --package-manager          : Select package manager for dependencies [ find | cmake | find-or-cmake | conan ] (default = find)
-f | --extra-flags [=arg]       : Extra CMake flags (defailt = none)
-g | --compiler [=arg]          : Compiler        | GNU | Clang | Tau (default = "")
-G | --generator [=arg]         : CMake generator  | many options... | (default = "CodeBlocks - Unix Makefiles")
   | --gcc-toolchain [=arg]     : Path to GCC toolchain. Use with Clang if it can't find stdlib (defailt = none)
-h | --help                     : Help. Shows this text.
-i | --install-prefix   [=path] : Install directory of DMRG (default = CMAKE_INSTALL_PREFIX)
   | --install-pkg-dir [=path] : Install directory of dependencies (default = CMAKE_INSTALL_PREFIX)
-j | --make-threads [=num]      : Number of threads used by Make build (default = 8)
-l | --clear-libs [=args]       : Clear libraries in comma separated list 'lib1,lib2...'. "all" deletes all.
-s | --enable-shared            : Enable shared library linking (default is static)
     --shared [=ON/OFF]         : Alternative to --enable-shared (default is OFF).
   | --enable-threads           : Enable C++11 threading used in Eigen::Tensor
   | --enable-mkl               : Enable Intel MKL
   | --enable-lto               : Enable Link Time Optimization
   | --enable-asan              : Enable runtime sanitizers, i.e. -fsanitize=address
   | --enable-pch               : Enable Precompiled Headers
   | --enable-ccache            : Enable ccache for speeding up compilation
   | --enable-coverage          : Enable test coverage
-t | --target [=args]           : Select build target [ CMakeTemplate | all-tests | test-<name> ]  (default = none)
   | --enable-tests             : Enable CTest (builds test targets)
   | --run-tests                : Run CTest
   | --print-cmake-error        : Prints CMakeError.log upon failure
   | --no-modules               : Disable use of "module load"
-v | --verbose                  : Verbose CMake and Makefiles
   | --verbose-make             : Verbose Makefiles
   | --verbose-cmake            : Verbose CMake
EXAMPLE:
./build.sh --arch native -b Release  --make-threads 8   --enable-shared --package-manager=find
EOF
  exit 1
}


# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"   -o ha:b:cl:df:g:G:j:s:t:v \
                --long "\
                help\
                arch:\
                build-type:\
                target:\
                clear-cmake\
                clear-libs:\
                compiler:\
                dry-run\
                default-tetralith\
                default-kraken\
                install-prefix:\
                install-pkg-dir:\
                package-manager:\
                enable-shared\
                shared:\
                gcc-toolchain:\
                make-threads:\
                enable-threads\
                enable-mkl\
                enable-lto\
                enable-asan\
                enable-coverage\
                enable-pch\
                enable-ccache\
                enable-tests\
                run-tests\
                no-modules\
                print-cmake-error
                verbose\
                verbose-make\
                verbose-cmake\
                generator:\
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
package_manager="find"
enable_threads="OFF"
enable_mkl="OFF"
enable_lto="OFF"
enable_asan="OFF"
enable_coverage="OFF"
enable_ccache="OFF"
enable_pch="OFF"
enable_tests="OFF"
run_tests="OFF"
install_prefix="~/dmrg-install"
install_pkg_dir="~/dmrg-deps-install"
make_threads=8
print_cmake_error="OFF"
verbose="OFF"
verbose_make="OFF"
verbose_cmake="OFF"
generator="Ninja"
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
       --default-tetralith)
                                    march="native"                  ; echo " * Architecture             : $march"          ;
                                    build_type="Release"            ; echo " * Build type               : $build_type"     ;
                                    compiler="GCC"                  ; echo " * C++ Compiler             : $compiler"       ;
                                    package_manager="conan"         ; echo " * Package Manager          : $package_manager";
                                    enable_shared="OFF"             ; echo " * Shared libraries         : $enable_shared"  ;
                                    enable_threads="ON"             ; echo " * C++11 Threads            : $enable_threads" ;
                                    enable_mkl="ON"                 ; echo " * Intel MKL                : $enable_mkl"     ;
                                    enable_lto="ON"                 ; echo " * Link Time Optimization   : $enable_lto"     ;
                                    enable_pch="ON"                 ; echo " * Preompiled headers       : $enable_pch"     ;
                                    enable_ccache="ON"              ; echo " * Ccache                   : $enable_ccache"  ;
                                    make_threads=16                 ; echo " * MAKE threads             : $make_threads"   ;
                                    target="all"                    ; echo " * CMake Build target       : $target"         ;
                                    verbose_cmake="ON"              ; echo " * Verbose cmake            : $verbose_cmake"  ;
                                    shift ;;
       --default-kraken)
                                    march="haswell"                 ; echo " * Architecture             : $march"          ;
                                    build_type="Release"            ; echo " * Build type               : $build_type"     ;
                                    compiler="GCC"                  ; echo " * C++ Compiler             : $compiler"       ;
                                    package_manager="conan"         ; echo " * Package Manager          : $package_manager";
                                    enable_shared="OFF"             ; echo " * Shared libraries         : $enable_shared"  ;
                                    enable_threads="ON"             ; echo " * C++11 Threads            : $enable_threads" ;
                                    enable_mkl="ON"                 ; echo " * Intel MKL                : $enable_mkl"     ;
                                    enable_lto="ON"                 ; echo " * Link Time Optimization   : $enable_lto"     ;
                                    enable_pch="ON"                 ; echo " * Preompiled headers       : $enable_pch"     ;
                                    enable_ccache="ON"              ; echo " * Ccache                   : $enable_ccache"  ;
                                    make_threads=32                 ; echo " * MAKE threads             : $make_threads"   ;
                                    target="all"                    ; echo " * CMake Build target       : $target"         ;
                                    verbose_cmake="ON"              ; echo " * Verbose cmake            : $verbose_cmake"  ;
                                    shift ;;
    -i|--install-prefix)            install_prefix=$2               ; echo " * Install Prefix           : $2"      ; shift 2 ;;
       --install-pkg-dir)           install_pkg_dir=$2              ; echo " * Install Pkg dir          : $2"      ; shift 2 ;;
       --package-manager)           package_manager=$2              ; echo " * Package Manager          : $2"      ; shift 2 ;;
    -f|--extra-flags)               extra_flags=$2                  ; echo " * Extra CMake flags        : $2"      ; shift 2 ;;
    -g|--compiler)                  compiler=$2                     ; echo " * C++ Compiler             : $2"      ; shift 2 ;;
    -G|--generator)                 generator=$2                    ; echo " * CMake generator          : $2"      ; shift 2 ;;
    -j|--make-threads)              make_threads=$2                 ; echo " * MAKE threads             : $2"      ; shift 2 ;;
    -s|--enable-shared)             enable_shared="ON"              ; echo " * Shared libraries         : ON"      ; shift   ;;
       --shared)                    enable_shared=$2                ; echo " * Shared libraries         : $2"      ; shift 2 ;;
    -t|--target)                    target=$2                       ; echo " * CMake Build target       : $2"      ; shift 2 ;;
       --enable-threads)            enable_threads="ON"             ; echo " * C++11 Threads            : ON"      ; shift   ;;
       --enable-mkl)                enable_mkl="ON"                 ; echo " * Intel MKL                : ON"      ; shift   ;;
       --enable-lto)                enable_lto="ON"                 ; echo " * Link Time Optimization   : ON"      ; shift   ;;
       --enable-asan)               enable_asan="ON"                ; echo " * Runtime sanitizers       : ON"      ; shift   ;;
       --enable-coverage)           enable_coverage="ON"            ; echo " * Coverage                 : ON"      ; shift   ;;
       --enable-pch)                enable_pch="ON"                 ; echo " * Preompiled headers       : ON"      ; shift   ;;
       --enable-ccache)             enable_ccache="ON"              ; echo " * Ccache                   : ON"      ; shift   ;;
       --enable-tests)              enable_tests="ON"               ; echo " * Build Tests              : ON"      ; shift   ;;
       --run-tests)                 run_tests="ON"                  ; echo " * Run CTest                : ON"      ; shift   ;;
       --no-modules)                no_modules="ON"                 ; echo " * Disable module load      : ON"      ; shift   ;;
       --print-cmake-error)         print_cmake_error="ON"          ; echo " * Print CMakeError.log     : ON"      ; shift   ;;
    -v|--verbose)                   verbose="ON"                    ; echo " * Verbose cmake & make     : ON"      ; shift   ;;
       --verbose-make)              verbose_make="ON"               ; echo " * Verbose make             : ON"      ; shift   ;;
       --verbose-cmake)             verbose_cmake="ON"              ; echo " * Verbose cmake            : ON"      ; shift   ;;
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

if [[ "$verbose" =~ ON|on|On|TRUE|True|true ]]; then
    verbose_make="ON"
    verbose_cmake="ON"
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

if [[ ! "$package_manager" =~ find|cmake|conan ]]; then
    echo "Package manager unsupported: $package_manager"
    exit 1
fi


if [ -n "$CONDA_PREFIX" ] ; then
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
        module load foss/2020b
        module load CMake
        if [ "$enable_mkl" = "ON" ] ; then
            export MKLROOT=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a/mkl
            export EBROOTIMKL=/software/sse/easybuild/prefix/software/imkl/2019.1.144-iimpi-2019a
        elif [[ "$package_manager" =~ find ]]; then
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
        if [[ "$package_manager" =~ find ]] ; then
                module load HDF5
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
echo " * CC                       :   $CC $($CC -dumpversion) at $(which $CC)"
echo " * CXX                      :   $CXX $($CXX -dumpversion) at $(which $CXX)"
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
          -DBUILD_SHARED_LIBS:BOOL=$enable_shared
          -DCMAKE_INSTALL_PREFIX:PATH=$install_prefix
          -DDMRG_DEPS_INSTALL_DIR:PATH=$install_pkg_dir
          -DCMAKE_VERBOSE_MAKEFILE=$verbose_make
          -DDMRG_PRINT_INFO=$verbose_cmake
          -DDMRG_PRINT_CHECKS=$verbose_cmake
          -DDMRG_PACKAGE_MANAGER=$package_manager
          -DDMRG_MICROARCH=$march
          -DDMRG_ENABLE_TESTS:BOOL=$enable_tests
          -DDMRG_ENABLE_THREADS=$enable_threads
          -DDMRG_ENABLE_MKL=$enable_mkl
          -DDMRG_ENABLE_LTO=$enable_lto
          -DDMRG_ENABLE_ASAN=$enable_asan
          -DDMRG_ENABLE_COVERAGE=$enable_coverage
          -DDMRG_ENABLE_PCH:BOOL=$enable_pch
          -DDMRG_ENABLE_CCACHE:BOOL=$enable_ccache
          $extra_flags
           -G $generator
           ../../
    cmake --build . --target $target --parallel $make_threads
EOF

if [ -z "$dry_run" ] ;then
    cmake -E make_directory build/$build_type
    cd build/$build_type
    cmake -DCMAKE_BUILD_TYPE=$build_type \
          -DBUILD_SHARED_LIBS:BOOL=$enable_shared \
          -DCMAKE_INSTALL_PREFIX:PATH=$install_prefix \
          -DDMRG_DEPS_INSTALL_DIR:PATH=$install_pkg_dir \
          -DCMAKE_VERBOSE_MAKEFILE=$verbose_make \
          -DDMRG_PRINT_INFO=$verbose_cmake \
          -DDMRG_PRINT_CHECKS=$verbose_cmake \
          -DDMRG_PACKAGE_MANAGER=$package_manager \
          -DDMRG_MICROARCH=$march \
          -DDMRG_ENABLE_TESTS:BOOL=$enable_tests \
          -DDMRG_ENABLE_THREADS=$enable_threads \
          -DDMRG_ENABLE_MKL=$enable_mkl \
          -DDMRG_ENABLE_LTO=$enable_lto \
          -DDMRG_ENABLE_ASAN=$enable_asan \
          -DDMRG_ENABLE_COVERAGE=$enable_coverage \
          -DDMRG_ENABLE_PCH:BOOL=$enable_pch \
          -DDMRG_ENABLE_CCACHE:BOOL=$enable_ccache \
          $extra_flags \
           -G "$generator" \
           ../../
    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo ""
            echo "Exit code: $exit_code"
            if [ "$print_cmake_error" = "ON" ] ; then
                echo "==================================================================="
                echo "CMakeFiles/CMakeError.log:"
                cat CMakeFiles/CMakeError.log
                echo "==================================================================="
            fi
            exit "$exit_code"
    fi

    cmake --build . --target $target --parallel $make_threads
    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo ""
            echo "Exit code: $exit_code"
            if [ "$print_cmake_error" = "ON" ] ; then
                echo "==================================================================="
                echo "CMakeFiles/CMakeError.log:"
                cat CMakeFiles/CMakeError.log
                echo "==================================================================="
            fi
            exit "$exit_code"
    fi


    if [ "$run_tests" = "ON" ] ;then
        if [[ "$target" == *"test-"* ]]; then
            ctest --build-config $build_type --verbose  --output-on-failure --build-target $target -R $target
        else
            ctest --build-config $build_type --verbose  --output-on-failure -R dmrg
        fi
    fi

    exit_code=$?
    if [ "$exit_code" != "0" ]; then
            echo "Exit code: $exit_code"
            exit "$exit_code"
    fi

fi


