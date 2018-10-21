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
enable_language(Fortran)

include(ExternalProject)
ExternalProject_Add(LAPACK
    GIT_REPOSITORY      https://github.com/Reference-LAPACK/lapack.git
    GIT_TAG             master
    GIT_PROGRESS 1
    PREFIX              "${INSTALL_DIRECTORY}/lapack"

	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	-DBLA_VENDOR=${BLA_VENDOR}
	-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
	-DCMAKE_C_FLAGS=-w
	-DBUILD_SHARED_LIBS=ON #Otherwise it wont work with the new arpack-ng
	-DUSE_OPTIMIZED_BLAS:bool=ON
	)
ExternalProject_Get_Property(LAPACK install_dir)

#ExternalProject_Get_Property(LAPACK install_dir)




