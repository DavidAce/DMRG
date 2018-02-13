#This script is supposed to be run from an execute_script() command such as:
#
#
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -E make_directory build
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/arpackpp)
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND}
#                    -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                    -G ${CMAKE_GENERATOR} ../
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/arpackpp/build )
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} --build .
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/arpackpp/build )
#
#            include(${PROJECT_SOURCE_DIR}/libs/arpackpp/FindArpackpp.cmake)
#            include_directories(${ARPACKPP_INCLUDE_DIR})
#
#
# It will generate a file called FindEigen3.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/arpackpp/FindArpackpp.cmake
#



cmake_minimum_required(VERSION 3.7)

project(ARPACKPP-download NONE)
include(ExternalProject)
ExternalProject_Add(ARPACKPP
        GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
        GIT_TAG             master
        PREFIX              "${INSTALL_DIRECTORY}/arpackpp"
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
        -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
        -DARPACK_LIB:PATH=${ARPACK_LIBRARIES}
        -DCMAKE_C_FLAGS=-w
        )



