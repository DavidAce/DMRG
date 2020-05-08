cmake_minimum_required(VERSION 3.12)
#include(../GetNumThreads.cmake)
#get_num_threads(NUM_THREADS)
#set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${NUM_THREADS})

project(external-glog)

include(ExternalProject)
ExternalProject_Add(external_glog
        URL https://github.com/google/glog/archive/v0.4.0.tar.gz
        URL_MD5 "0daea8785e6df922d7887755c3d100d0"
        BUILD_ALWAYS TRUE
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -Dgflags_DIR:BOOL:PATH=${gflags_DIR}
        -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )