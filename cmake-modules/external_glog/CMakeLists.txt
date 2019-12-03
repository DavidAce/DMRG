cmake_minimum_required(VERSION 3.10)
project(external-glog)

include(ExternalProject)
ExternalProject_Add(external_glog
        GIT_REPOSITORY https://github.com/google/glog.git
        GIT_TAG v0.4.0
        GIT_PROGRESS false
        GIT_SHALLOW true
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}/glog
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_FLAGS:STRING=${extra_flags}
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DCMAKE_FIND_LIBRARY_SUFFIXES=${CMAKE_STATIC_LIBRARY_SUFFIX}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )