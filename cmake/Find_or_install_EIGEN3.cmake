# The following script will find a native Eigen3 installation on your system.
#
# If these packages are not found, they will be downloaded from a repository and installed
# Either way, this script will define the following variables:
#

# EIGEN3_FOUND - system has eigen lib with correct version
# EIGEN3_INCLUDE_DIR - the eigen include directory
# EIGEN3_VERSION - eigen version
#
# After execution, it is enough to include directory:
#   include_directory(${EIGEN3_INCLUDE_DIR}), which is done here.
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_EIGEN3.cmake)


message("SEARCHING FOR PRE-INSTALLED LIBRARIES: EIGEN3")
find_package(Eigen3 3.3)                                ### Find and define includes for Eigen Library
if(EIGEN3_FOUND)
    message("FOUND PRE-INSTALLED EIGEN3:   ${EIGEN3_INCLUDE_DIR}")
else()
    message("DOWNLOADING EIGEN3...")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/eigen3
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
                -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
                -G ${CMAKE_GENERATOR} ../../eigen3
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3)
    execute_process(
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "Eigen3 not found and failed to install: ${res_var}")
    endif()
    include(${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake)
    include_directories(${EIGEN3_INCLUDE_DIR})
    message("SUCCESSFULLY INSTALLED EIGEN3:   ${EIGEN3_INCLUDE_DIR}")
    message("BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3/log_build.txt")
endif()


