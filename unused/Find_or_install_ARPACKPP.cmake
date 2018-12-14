# The following script will find a native Arpackpp installation on your system.
#
# If these packages are not found, they will be downloaded from a repository and installed
# Either way, this script will define the following variables:
#

# ARPACKPP_FOUND - system has eigen lib with correct version
# ARPACKPP_INCLUDE_DIR - the eigen include directory
# ARPACKPP_VERSION - eigen version
#
# After execution, it is enough to include directory:
#   include_directory(${ARPACKPP_INCLUDE_DIR}), which is done here.
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_ARPACKPP.cmake)


message("SEARCHING FOR LIBRARY: ARPACKPP")
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/arpackpp/FindArpackpp.cmake")
    include(${PROJECT_SOURCE_DIR}/libs/arpackpp/FindArpackpp.cmake)
    message("FOUND PREVIOUSLY INSTALLED ARPACKPP:   ${ARPACKPP_INCLUDE_DIR}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACKPP_INCLUDE_DIR})
    return()
endif()

message(STATUS "SEARCHING FOR ARPACKPP IN SYSTEM...")
#find_library(ARPACKPP_LIBRARIES
#        NAMES "arpack++" "arpackpp"
#        PATH_SUFFIXES "lib" "lib32" "lib64"
#        )
#find_path(ARPACKPP_INCLUDE_DIR
#        NAMES "arssym.h"
#        PATH_SUFFIXES "include" "include/arpack++"
#        )
#find_library(ARPACKPP_INCLUDE_DIR NAMES libarpack++ arpackpp libarpack++2-dev libarpack++2c2a )                                ### Find and define includes for Eigen Library
if(ARPACKPP_LIBRARIES AND ARPACKPP_INCLUDE_DIR)
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACKPP_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${ARPACKPP_LIBRARIES})
    get_cmake_property(_variableNames VARIABLES)
    return()
else()
    message(STATUS "DOWNLOADING ARPACKPP...")
    set(INSTALL_DIRECTORY ../libs)
    if(BLAS_LIBRARIES AND LAPACK_LIBRARIES AND ARPACK_LIBRARIES)
        execute_process(
                COMMAND ${CMAKE_COMMAND} -E make_directory tmp/arpackpp
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
        execute_process(
                COMMAND ${CMAKE_COMMAND}
                    -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
                    -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
                    -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
                    -DARPACK_LIBRARIES:PATH=${ARPACK_LIBRARIES}
                    -G ${CMAKE_GENERATOR} ../../arpackpp
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp)
        execute_process(
                OUTPUT_QUIET
                OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp/log_build.txt
                COMMAND ${CMAKE_COMMAND} --build .
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp
                RESULT_VARIABLE res_var)
        if(NOT "${res_var}" STREQUAL "0")
            message(FATAL_ERROR "ARPACKPP not found and failed to install: ${res_var}")
        endif()
        # Generate a file that can be detected if the library was installed successfully
        find_path(ARPACKPP_INCLUDE_DIR
                HINTS ${INSTALL_DIRECTORY}/arpackpp
                NO_DEFAULT_PATH NO_MODULE
                NAMES "arssym.h"
                PATH_SUFFIXES "include" "include/arpack++" "include/arpackpp"
                REQUIRED
                )
#        set(ARPACKPP_INC_DIR ${INSTALL_DIRECTORY}/arpackpp/include ${INSTALL_DIRECTORY}arpackpp/include/arpackpp)
        file(WRITE ${INSTALL_DIRECTORY}/arpackpp/FindArpackpp.cmake "set(ARPACKPP_INCLUDE_DIR   ${ARPACKPP_INCLUDE_DIR})\n")
        # Include that file
        include(${INSTALL_DIRECTORY}/arpackpp/FindArpackpp.cmake)
        target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACKPP_INCLUDE_DIR})
        message(STATUS "SUCCESSFULLY INSTALLED ARPACKPP:   ${ARPACKPP_INCLUDE_DIR}")
        message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp/log_build.txt")
    else()
        message(FATAL_ERROR "Please install LAPACK, BLAS and ARPACK before ARPACKPP.")

    endif()
endif()


