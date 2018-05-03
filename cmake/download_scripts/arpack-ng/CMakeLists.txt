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
#set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "-Wl --no-as-needed")
include(ExternalProject)
ExternalProject_Add(ARPACK
	GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
#	GIT_TAG             3.5.0
	GIT_TAG             master
	PREFIX              "${INSTALL_DIRECTORY}/arpack-ng"
	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
	-DCMAKE_C_FLAGS=-w
	-DEXAMPLES=OFF
	-DCOVERALLS=OFF
	-DCMAKE_BUILD_TYPE=Release
	-DMPI=OFF
	-DBUILD_SHARED_LIBS=OFF
	-DICB=OFF
	-DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
	-DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
	-DEXTRA_LDLAGS=mkl_gf_lp64
	)

