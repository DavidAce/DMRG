cmake_minimum_required(VERSION 3.10)
project(external-h5pp)
include(FetchContent)
message("EXTERNAL_BUILD_DIR: ${EXTERNAL_BUILD_DIR}")
message("EXTERNAL_INSTALL_DIR: ${EXTERNAL_INSTALL_DIR}")

FetchContent_Declare(h5pp
        GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
        GIT_TAG master
        GIT_PROGRESS 1
        GIT_SHALLOW 1
        PREFIX      ${EXTERNAL_BUILD_DIR}/h5pp
        INSTALL_DIR ${EXTERNAL_INSTALL_DIR}/h5pp
#        BUILD_ALWAYS OFF
#        UPDATE_COMMAND ""
#        TEST_COMMAND ${CMAKE_MAKE_PROGRAM} test
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_STANDARD=17
        -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=ON
        -DCMAKE_CXX_EXTENSIONS:BOOL=OFF
        -DENABLE_TESTS:BOOL=ON
        -DTHIRD_PARTY_DIR:PATH=<INSTALL_DIR>
        -DDOWNLOAD_MISSING:BOOL=ON
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console, like up-to-date and installing
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )

FetchContent_GetProperties(h5pp)
if(NOT h5pp_POPULATED)
    message(STATUS "Fetching content")
    FetchContent_Populate(h5pp)
    add_subdirectory(${h5pp_SOURCE_DIR} ${h5pp_BINARY_DIR})
    message(STATUS "h5pp_SOURCE_DIR ${h5pp_SOURCE_DIR}")
    message(STATUS "h5pp_BINARY_DIR ${h5pp_BINARY_DIR}")

endif()

include(cmake-modules/PrintTargetInfo.cmake)
print_target_info(h5pp)
print_target_info(h5pp::h5pp)
