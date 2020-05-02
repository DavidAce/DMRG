cmake_minimum_required(VERSION 3.12)
#include(../GetNumThreads.cmake)
#get_num_threads(NUM_THREADS)
#set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${NUM_THREADS})

project(external-glog)

include(ExternalProject)
ExternalProject_Add(external_glog
        GIT_REPOSITORY https://github.com/google/glog.git
        GIT_TAG v0.4.0
        GIT_PROGRESS false
        GIT_SHALLOW true
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
#        BUILD_COMMAND ${CMAKE_COMMAND} --build . -- -j ${NUM_THREADS}
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#        -DCMAKE_CXX_STANDARD=11
        -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=ON
        -DCMAKE_CXX_EXTENSIONS:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -Dgflags_DIR:BOOL:PATH=${gflags_DIR}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )