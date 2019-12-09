cmake_minimum_required(VERSION 3.10)
project(external-gflags)

include(ExternalProject)
ExternalProject_Add(external_gflags
        GIT_REPOSITORY https://github.com/gflags/gflags.git
        GIT_TAG v2.2.2
        GIT_PROGRESS false
        GIT_SHALLOW true
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}/gflags
        UPDATE_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DGFLAGS_BUILD_SHARED_LIBS:BOOL=OFF
        -DGFLAGS_BUILD_STATIC_LIBS:BOOL=ON
#        -DBUILD_gflags_nothreads_LIBS:BOOL=OFF
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        )