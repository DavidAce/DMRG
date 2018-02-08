#This script is supposed to be run from an execute_script() command such as:
#
#
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -E make_directory build
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND}
#                                   -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                                   -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
#                                   -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
#                                   -G ${CMAKE_GENERATOR} ../
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
# Note that you must pass the two parameters so that arpack can find these two dependencies
#   	-DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
#       -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
# It will generate a file called FindArpack.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/arpack/FindArpack.cmake
#



cmake_minimum_required(VERSION 3.7)

project(ARPACK-download NONE)
include(ExternalProject)
ExternalProject_Add(ARPACK
	GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
	GIT_TAG             master
	PREFIX              "${INSTALL_DIRECTORY}/arpack-ng"
	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
	-DCMAKE_C_FLAGS=-w
	-DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
	-DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
        )

ExternalProject_Get_Property(ARPACK install_dir)
set(ARPACK_LIB_DIR ${install_dir}/lib)
set(ARPACK_INC_DIR ${install_dir}/include)
set(ARPACK_BIN_DIR ${install_dir}/bin)
execute_process(
		COMMAND  find ${install_dir} -name arpack-ng-config.cmake -exec dirname {} \; -quit
		OUTPUT_VARIABLE ARPACK_CMAKE_DIR
)
message("ARPACK_CMAKE_DIR: ${ARPACK_CMAKE_DIR}")
file(WRITE  ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_LIB_DIR   ${ARPACK_LIB_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_INC_DIR   ${ARPACK_INC_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_BIN_DIR   ${ARPACK_BIN_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_CMAKE_DIR   ${ARPACK_CMAKE_DIR})\n")





