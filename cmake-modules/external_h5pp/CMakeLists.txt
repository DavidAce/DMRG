cmake_minimum_required(VERSION 3.10)
project(external-h5pp)

include(ExternalProject)
ExternalProject_Add(external_h5pp
        GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
        GIT_TAG v0.4.3
        GIT_PROGRESS 1
        GIT_SHALLOW 1
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        TEST_COMMAND ${CMAKE_MAKE_PROGRAM} test
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_STANDARD=17
        -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=ON
        -DCMAKE_CXX_EXTENSIONS:BOOL=OFF
        -DENABLE_TESTS:BOOL=ON
        -DDOWNLOAD_MISSING:BOOL=ON
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DAPPEND_LIBSUFFIX:BOOL=ON
        -DH5PP_DIRECTORY_HINTS:PATH=${H5PP_DIRECTORY_HINTS}
        -Dhdf5_DIR:PATH=${hdf5_DIR}
        -DEigen3_DIR:PATH=${Eigen3_DIR}
        -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console, like up-to-date and installing
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )