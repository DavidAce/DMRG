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



find_library(ARPACKPP_INCLUDE_DIR NAMES libarpack++ arpackpp libarpack++2-dev libarpack++2c2a )                                ### Find and define includes for Eigen Library
if(ARPACKPP_INCLUDE_DIR)
    message("FOUND NATIVE ARPACKPP:   ${ARPACKPP_INCLUDE_DIR}")
else()
    message("DOWNLOADING ARPACKPP...")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/arpackpp
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
                -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
                -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
                -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
                -DARPACK_LIBRARIES:PATH=${ARPACK_LIBRARIES}
                -G ${CMAKE_GENERATOR} ../../arpackpp
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp)
    execute_process(
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpackpp
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "ARPACKPP not found and failed to install: ${res_var}")
    endif()
    include(${PROJECT_SOURCE_DIR}/libs/arpackpp/FindArpackpp.cmake)
    include_directories(${ARPACKPP_INCLUDE_DIR})
    message("SUCCESSFULLY INSTALLED ARPACKPP:   ${ARPACKPP_INCLUDE_DIR}")
endif()


