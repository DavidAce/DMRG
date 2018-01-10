#This script is supposed to be run from an execute_script() command such as:
#
#
#        execute_process(
#                COMMAND ${CMAKE_COMMAND} -E make_directory build
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/eigen3)
#        execute_process(
#                COMMAND ${CMAKE_COMMAND}
#                -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                -G ${CMAKE_GENERATOR} ../
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/eigen3/build )
#        execute_process(
#                COMMAND ${CMAKE_COMMAND} --build .
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/eigen3/build )
#        
#        include(${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake)
#        include_directories(${EIGEN3_INCLUDE_DIR})
#        message("SUCCESSFULLY INSTALLED EIGEN3:   ${EIGEN3_INCLUDE_DIR}")
#
#
# It will generate a file called FindEigen3.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake
#

#        URL "https://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2"
#        GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
#        HG_REPOSITORY https://bitbucket.org/eigen/eigen


cmake_minimum_required(VERSION 3.7)
project(EIGEN3-download NONE)
include(ExternalProject)
ExternalProject_Add(EIGEN3
        GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX "${INSTALL_DIRECTORY}/eigen3"
        TEST_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        -Wno-dev
        )

ExternalProject_Get_Property(EIGEN3 install_dir)
set(EIGEN3_INC_DIR ${install_dir}/include/eigen3)
file(WRITE ${INSTALL_DIRECTORY}/eigen3/FindEigen3.cmake "set(EIGEN3_INCLUDE_DIR   ${EIGEN3_INC_DIR})\n")

