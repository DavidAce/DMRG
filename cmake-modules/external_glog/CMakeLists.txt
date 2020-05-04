cmake_minimum_required(VERSION 3.12)
#include(../GetNumThreads.cmake)
#get_num_threads(NUM_THREADS)
#set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${NUM_THREADS})

project(external-glog)

include(ExternalProject)
ExternalProject_Add(external_glog
        GIT_REPOSITORY https://github.com/google/glog.git
        GIT_TAG v0.4.0
        GIT_PROGRESS FALSE
        GIT_SHALLOW TRUE
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