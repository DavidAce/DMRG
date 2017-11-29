# The following script will find a native BLAS and LAPACK installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# ARPACK_LIBRARIES          - Path to liblapack.a
#
# ARPACK_LIB                - Path do lapack library directory
# ARPACK_INC                - Path to lapack include directory
# ARPACK_BIN                - Path to lapack binary directory
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${ARPACK_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_ARPACK.cmake)

find_library(ARPACK_LIBRARIES NAMES arpack arpack-ng libarpack libarpack2 libarpack2-dev)
if(ARPACK_LIBRARIES)
    message("FOUND NATIVE ARPACK:   ${ARPACK_LIBRARIES}")
elseif(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
    message("DOWNLOADING ARPACK...")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/arpack-ng
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
                -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
                -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
                -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
                -G ${CMAKE_GENERATOR} ../../arpack-ng
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng )
    execute_process(
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "ARPACK not found and failed to install: ${res_var}")
    endif()

    include(${PROJECT_SOURCE_DIR}/libs/arpack-ng/FindArpack.cmake)
    set(ARPACK_LIBRARIES "")
    get_libraries(${ARPACK_LIB_DIR} arpack  ARPACK_LIBRARIES)
    include_directories(${ARPACK_INC_DIR})
    message("SUCCESSFULLY INSTALLED ARPACK:   ${ARPACK_LIBRARIES}")
else()

    message("Please install LAPACK and BLAS before ARPACK.")

endif()


