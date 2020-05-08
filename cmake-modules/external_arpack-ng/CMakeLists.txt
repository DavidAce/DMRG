cmake_minimum_required(VERSION 3.12)
project(external-arpack-ng)

include(ExternalProject)
ExternalProject_Add(external_ARPACK-NG
        URL https://github.com/opencollab/arpack-ng/archive/3.7.0.tar.gz
        URL_MD5 "6fc6c6bf78dbd4f144595ef0675c8430"
        BUILD_ALWAYS TRUE
        PREFIX           ${CMAKE_BINARY_DIR}
        INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}
        UPDATE_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        #-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        -DCMAKE_C_FLAGS=${ARPACK_FLAGS}
        -DCMAKE_Fortran_FLAGS=${ARPACK_FLAGS}
        -DEXAMPLES=OFF
        -DMPI=OFF
        -DINTERFACE64=OFF
        -DBLAS_LIBRARIES=${BLAS_LIBRARIES}
        -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
        -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
        )