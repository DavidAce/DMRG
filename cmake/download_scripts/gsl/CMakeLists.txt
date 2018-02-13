#This script is supposed to be run from an execute_script() command such as:
#

#        find_package(GSL)
#        if (GSL_FOUND)
#            include_directories(${GSL_INCLUDE_DIRS})
#            message("FOUND NATIVE GSL:   ${GSL_LIBRARIES}")
#        else()
#            message("DOWNLOADING GSL...")
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -E make_directory build
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/gsl)
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND}
#                    -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                    -G ${CMAKE_GENERATOR} ../
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/gsl/build )
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} --build .
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/gsl/build )
#            include(${PROJECT_SOURCE_DIR}/libs/gsl/FindGSL.cmake)
#            include_directories(${GSL_INC_DIR})
#            message("SUCCESSFULLY INSTALLED GSL:   ${GSL_LIBRARIES}")
#        endif()
# It will generate a file called FindGSL.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/gsl/FindGSL.cmake
#



cmake_minimum_required(VERSION 3.7)
project(GSL-download NONE)
include(ExternalProject)
ExternalProject_Add(GSL
        URL      http://ftp.acc.umu.se/mirror/gnu.org/gnu/gsl/gsl-2.4.tar.gz
#        GIT_TAG             master
#        GIT_PROGRESS 1
#        UPDATE_DISCONNECTED 1
        PREFIX              "${INSTALL_DIRECTORY}/gsl"
        CONFIGURE_COMMAND
            cd <SOURCE_DIR> &&
            pwd &&
            ./configure CFLAGS=-w --quiet --silent --enable-silent-rules --enable-shared=no --prefix=<INSTALL_DIR>
        BUILD_COMMAND
            cd <SOURCE_DIR> &&
            pwd &&
            ${CMAKE_MAKE_PROGRAM} --quiet --silent
        INSTALL_COMMAND
            cd <SOURCE_DIR> &&
            pwd &&
            ${CMAKE_MAKE_PROGRAM} --quiet install
        )






