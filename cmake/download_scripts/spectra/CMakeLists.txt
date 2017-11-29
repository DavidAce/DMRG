#This script is supposed to be run from an execute_script() command such as:
#
#
#        execute_process(
#                COMMAND ${CMAKE_COMMAND} -E make_directory build
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/spectra)
#        execute_process(
#                COMMAND ${CMAKE_COMMAND}
#                -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                -G ${CMAKE_GENERATOR} ../
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/spectra/build )
#        execute_process(
#                COMMAND ${CMAKE_COMMAND} --build .
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/spectra/build )
#        
#        include(${PROJECT_SOURCE_DIR}/libs/spectra/FindSpectra.cmake)
#        include_directories(${SPECTRA_INCLUDE_DIR})
#        message("SUCCESSFULLY INSTALLED SPECTRA:   ${SPECTRA_INCLUDE_DIR}")
#
#
# It will generate a file called FindSpectra.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/spectra/FindSpectra.cmake
#

#        URL "https://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2"
#        HG_REPOSITORY https://bitbucket.org/eigen/eigen


cmake_minimum_required(VERSION 3.8)
project(SPECTRA-download NONE)
include(ExternalProject)
ExternalProject_Add(SPECTRA
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX ${INSTALL_DIRECTORY}/spectra
        TEST_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        -Wno-dev
        )

ExternalProject_Get_Property(SPECTRA install_dir)
set(SPECTRA_INC_DIR ${install_dir}/src/SPECTRA/include)
file(WRITE ${INSTALL_DIRECTORY}/spectra/FindSpectra.cmake "set(SPECTRA_INCLUDE_DIR   ${SPECTRA_INC_DIR})\n")

