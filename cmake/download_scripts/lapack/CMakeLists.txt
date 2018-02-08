#This script is supposed to be run from an execute_script() command such as:
#
#
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -E make_directory build
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs -G ${CMAKE_GENERATOR} ../
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/build )
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} --build .
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/build)
#
#            include(${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake)
#            include_directories(${Lapack_INC_DIR})
#            message("SUCCESSFULLY INSTALLED LAPACK:   ${LAPACK_LIBRARIES}")
#            message("SUCCESSFULLY INSTALLED BLAS:     ${BLAS_LIBRARIES}")
#
#
# It will generate a file called FindLapack.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake
#



cmake_minimum_required(VERSION 3.7)

project(LAPACK-download NONE)
include(ExternalProject)
ExternalProject_Add(LAPACK
    GIT_REPOSITORY      https://github.com/Reference-LAPACK/lapack.git
    GIT_TAG             master
    GIT_PROGRESS 1
    PREFIX              "${INSTALL_DIRECTORY}/lapack"
	CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
		-DCMAKE_C_FLAGS=-w
        )

ExternalProject_Get_Property(LAPACK install_dir)
set(LAPACK_LIB_DIR ${install_dir}/lib)
set(LAPACK_INC_DIR ${install_dir}/include)
set(LAPACK_BIN_DIR ${install_dir}/bin)
set(LAPACK_CMAKE_DIR ${install_dir}/lib/cmake/lapack-3.8.0)


file(WRITE ${INSTALL_DIRECTORY}/lapack/FindLapack.cmake "set(LAPACK_LIB_DIR    ${LAPACK_LIB_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/lapack/FindLapack.cmake "set(LAPACK_INC_DIR   ${LAPACK_INC_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/lapack/FindLapack.cmake "set(LAPACK_BIN_DIR   ${LAPACK_BIN_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/lapack/FindLapack.cmake "set(LAPACK_CMAKE_DIR   ${LAPACK_CMAKE_DIR})\n")





