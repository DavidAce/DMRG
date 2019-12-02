cmake_minimum_required(VERSION 3.10)
project(external-h5pp)

include(ExternalProject)
ExternalProject_Add(external_h5pp
        GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX      ${EXTERNAL_BUILD_DIR}/h5pp
        INSTALL_DIR ${EXTERNAL_INSTALL_DIR}/h5pp
        BUILD_ALWAYS OFF
        UPDATE_COMMAND ""
        TEST_COMMAND ${CMAKE_MAKE_PROGRAM} test
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DENABLE_TESTS:BOOL=ON
        -DTHIRD_PARTY_DIR:PATH=${EXTERNAL_INSTALL_DIR}
        -DDOWNLOAD_MISSING:BOOL=ON
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console, like up-to-date and installing
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )