cmake_minimum_required(VERSION 3.10)
project(external-OpenBLAS)
include(../GetNumThreads.cmake)
get_num_threads(NUM_THREADS)
set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${NUM_THREADS})
message(STATUS "Building OpenBLAS with ${NUM_THREADS} threads")
include(ExternalProject)

ExternalProject_Add(external_OpenBLAS
        URL      https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz
        URL_MD5  5cd4ff3891b66a59e47af2d14cde4056
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/OpenBLAS
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        CONFIGURE_COMMAND ""
#
#        CMAKE_ARGS
#        -DTARGET:STRING=${TARGET}
#        -DDYNAMIC_ARCH:BOOL=OFF
#        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
#        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#        -DUSE_THREAD:BOOL=${USE_THREAD}
#        -DBINARY=64
#        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
#        -DCMAKE_C_FLAGS=-w
#        -DCMAKE_CXX_FLAGS=-w
#        -DCMAKE_FORTRAN_FLAGS=-w

#        BUILD_COMMAND export LD_LIBRARY_PATH=${GFORTRAN_PATH}:$LD_LIBRARY_PATH &&
#        export LDFLAGS="-L${GFORTRAN_PATH} $LDFLAGS" &&

        BUILD_IN_SOURCE 1
        BUILD_COMMAND
            export LD_LIBRARY_PATH=${GFORTRAN_PATH} &&
            export LDFLAGS=${LDFLAGS} &&
            export FFLAGS=${FFLAGS} &&
            ${CMAKE_MAKE_PROGRAM} -j ${NUM_THREADS}
                    TARGET=${TARGET}
                    DYNAMIC_ARCH=0
                    USE_THREAD=${USE_THREAD}
                    USE_OPENMP=${USE_OPENMP}
                    NO_AFFINITY=1
                    NO_WARMUP=1
                    QUIET_MAKE=0
                    NUM_THREADS=64
                    DEBUG=0
                    BINARY=64
                    GEMM_MULTITHREAD_THRESHOLD=32

        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM}  PREFIX=<INSTALL_DIR> install
        )