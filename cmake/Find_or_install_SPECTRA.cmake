# The following script will find a native Spectra installation on your system.
#
# If these packages are not found, they will be downloaded from a repository and installed
# Either way, this script will define the following variables:
#

# SPECTRA_FOUND - system has eigen lib with correct version
# SPECTRA_INCLUDE_DIR - the eigen include directory
# SPECTRA_VERSION - eigen version
#
# After execution, it is enough to include directory:
#   include_directory(${SPECTRA_INCLUDE_DIR}), which is done here.
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_SPECTRA.cmake)



message("DOWNLOADING SPECTRA...")
execute_process(
        COMMAND ${CMAKE_COMMAND} -E make_directory tmp/spectra
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
execute_process(
        COMMAND ${CMAKE_COMMAND}
        -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
        -G ${CMAKE_GENERATOR} ../../spectra
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/spectra)
execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/spectra
        RESULT_VARIABLE res_var)
if(NOT "${res_var}" STREQUAL "0")
    message(FATAL_ERROR "Spectra not found and failed to install: ${res_var}")
endif()
include(${PROJECT_SOURCE_DIR}/libs/spectra/FindSpectra.cmake)
include_directories(${SPECTRA_INCLUDE_DIR})
message("SUCCESSFULLY INSTALLED SPECTRA:   ${SPECTRA_INCLUDE_DIR}")



