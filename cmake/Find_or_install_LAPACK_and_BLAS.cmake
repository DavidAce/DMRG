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




message("SEARCHING FOR LIBRARIES: LAPACK and BLAS")
set(LAPACK_CMAKE_DIR_GUESS ${PROJECT_SOURCE_DIR}/libs/lapack/lib/cmake/lapack-3.8.0)
find_package(LAPACK PATHS "${LAPACK_CMAKE_DIR_GUESS}" NO_DEFAULT_PATH NO_MODULE)
if(LAPACK_FOUND)
    include(${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake)
    set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARIES})
    set(BLAS_LIBRARIES ${LAPACK_blas_LIBRARIES})
    include_directories(${LAPACK_INC_DIR})
    message(STATUS "FOUND PREVIOUSLY INSTALLED LAPACK:   ${LAPACK_LIBRARIES}")
    message(STATUS "FOUND PREVIOUSLY INSTALLED BLAS:     ${BLAS_LIBRARIES}")
    return()
endif()


message(STATUS "SEARCHING FOR LAPACK AND BLAS IN SYSTEM...")
find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "FOUND BLAS:     ${BLAS_LIBRARIES}")
    message(STATUS "FOUND LAPACK:   ${LAPACK_LIBRARIES}")
else()
    message(STATUS "DOWNLOADING BLAS AND LAPACK...")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E make_directory tmp/lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
	execute_process(
		COMMAND ${CMAKE_COMMAND} -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs -G ${CMAKE_GENERATOR} ../../lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack )
	execute_process(
        OUTPUT_QUIET
        OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack/log_build.txt
		COMMAND ${CMAKE_COMMAND} --build .
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack
        RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "LAPACK and BLAS not found and failed to install: ${res_var}")
    endif()
	include(${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake)
    set(LAPACK_LIBRARIES "")
    set(BLAS_LIBRARIES "")
    find_package(LAPACK PATHS "${LAPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
    set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARIES})
    set(BLAS_LIBRARIES ${LAPACK_blas_LIBRARIES})
    include_directories(${LAPACK_INC_DIR})
    message(STATUS "SUCCESSFULLY INSTALLED LAPACK:   ${LAPACK_LIBRARIES}")
    message(STATUS "SUCCESSFULLY INSTALLED BLAS:     ${BLAS_LIBRARIES}")
    message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack/log_build.txt")
endif()

#get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        if("${_variableName}" MATCHES "BLAS" OR "${_variableName}" MATCHES "LAPACK")
#            message(STATUS "${_variableName}=${${_variableName}}")
#        endif()
#    endforeach()