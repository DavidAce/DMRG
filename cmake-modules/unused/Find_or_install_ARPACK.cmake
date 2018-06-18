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

message("SEARCHING FOR LIBRARY: ARPACK")
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/arpack-ng/FindArpack.cmake" )
    include(${PROJECT_SOURCE_DIR}/libs/arpack-ng/FindArpack.cmake)
    find_package(arpack-ng PATHS "${ARPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
    set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES})
    set(ARPACK_INCLUDE_DIRS ${arpack_ng_INCLUDE_DIRS})
    if(arpack-ng_FOUND AND LAPACK_LIBRARIES AND BLAS_LIBRARIES)
        message(STATUS "FOUND PREVIOUSLY INSTALLED ARPACK:   ${ARPACK_LIBRARIES}")
        target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACK_INCLUDE_DIRS})
        target_link_libraries(${PROJECT_NAME} ${ARPACK_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${GFORTRAN_LIB})
        return()
    endif()
    message(WARNING "Found 'libs/arpack-ng/FindArpack.cmake' in , but failed to configure. You may need a clean build.")
endif()

message(STATUS "SEARCHING FOR ARPACK IN SYSTEM...")
find_library(ARPACK_LIBRARIES
        NAMES "arpack" "libarpack.a"
        PATH_SUFFIXES "lib" "lib32" "lib64"
)




if(ARPACK_LIBRARIES)
    message(STATUS "FOUND ARPACK:   ${ARPACK_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${ARPACK_LIBRARIES})
    return()
elseif(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
    message(STATUS "DOWNLOADING ARPACK...")
    message("Using BLAS_LIBRARIES  : ${BLAS_LIBRARIES}")
    message("Using LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
    set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)
    execute_process(COMMAND  export PATH="$PATH:/opt/intel/bin")
    execute_process(COMMAND  export LD_LIBRARY_PATH="$PATH:opt/intel/mkl/lib/intel64")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/arpack-ng
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng/log_commands.txt
            COMMAND
                ${CMAKE_COMMAND}
                -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
                -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
                -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
                -G ${CMAKE_GENERATOR} ../../arpack-ng
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng )
    execute_process(
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "ARPACK not found and failed to install: ${res_var}")
    endif()
    # Generate a file that can be detected if the library was installed successfully
    set(ARPACK_CMAKE_DIR ${INSTALL_DIRECTORY}/arpack-ng/lib/cmake)
    file(WRITE ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_CMAKE_DIR  ${ARPACK_CMAKE_DIR})\n")

    # Include that file
    include(${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake)
    find_package(arpack-ng PATHS "${ARPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
    set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES})
    set(ARPACK_INCLUDE_DIRS ${arpack_ng_INCLUDE_DIRS})

#









#    # Generate a file that can be detected if the library was installed successfully
#    set(ARPACK_CMAKE_DIR ${INSTALL_DIRECTORY}/arpack-ng)
#    file(WRITE ${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake "set(ARPACK_CMAKE_DIR  ${ARPACK_CMAKE_DIR})\n")
#
#    # Include that file
#    include(${INSTALL_DIRECTORY}/arpack-ng/FindArpack.cmake)
#    find_package(arpack-ng PATHS "${ARPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
#    set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES})
#    set(ARPACK_INCLUDE_DIRS ${arpack_ng_INCLUDE_DIRS})
#
#




    if(arpack-ng_FOUND)
        message(STATUS "SUCCESSFULLY INSTALLED ARPACK:   ${ARPACK_LIBRARIES}")
        message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng/log_build.txt")
        target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACK_INCLUDE_DIRS})
        target_link_libraries(${PROJECT_NAME} ${ARPACK_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${GFORTRAN_LIB})
    else()
    message(WARNING "ARPACK_LIBRARIES: ${ARPACK_LIBRARIES}")
    message(FATAL_ERROR "LAPACK or BLAS could not be found.")
    endif()
else()
    message(FATAL_ERROR "Please install LAPACK and BLAS before ARPACK.")
endif()

#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #        if("${_variableName}" MATCHES "Arpack")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #        endif()
#endforeach()