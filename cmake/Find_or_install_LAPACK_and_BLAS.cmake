# The following script will find a native BLAS and LAPACK installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# LAPACK_LIBRARIES          - Path to liblapack.a
# BLAS_LIBRARIES            - Path to libblas.a, which comes included in the lapack installation.
#
# LAPACK_LIB                - Path do lapack library directory
# LAPACK_INC                - Path to lapack include directory
# LAPACK_BIN                - Path to lapack binary directory
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_LAPACK_and_BLAS.cmake)

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    message("FOUND NATIVE BLAS:     ${BLAS_LIBRARIES}")
    message("FOUND NATIVE LAPACK:   ${LAPACK_LIBRARIES}")
else()
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E make_directory tmp/lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
	execute_process(
		COMMAND ${CMAKE_COMMAND} -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs -G ${CMAKE_GENERATOR} ../../lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack )
	execute_process(
		COMMAND ${CMAKE_COMMAND} --build .
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack
        RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "LAPACK and BLAS not found and failed to install: ${res_var}")
    endif()
	include(${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake)
    set(LAPACK_LIBRARIES "")
    set(BLAS_LIBRARIES "")
    get_libraries(${LAPACK_LIB_DIR} lapack  LAPACK_LIBRARIES)
    get_libraries(${LAPACK_LIB_DIR} blas    BLAS_LIBRARIES)
    include_directories(${LAPACK_INC_DIR})
    message("SUCCESSFULLY INSTALLED LAPACK:   ${LAPACK_LIBRARIES}")
    message("SUCCESSFULLY INSTALLED BLAS:     ${BLAS_LIBRARIES}")

endif()


