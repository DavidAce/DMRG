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
set(ARPACK_CMAKE_DIR_GUESS "${PROJECT_SOURCE_DIR}/libs/arpack-ng")
find_package(arpack-ng PATHS "${ARPACK_CMAKE_DIR_GUESS}" NO_DEFAULT_PATH NO_MODULE)
if(arpack-ng_FOUND)
    include(${PROJECT_SOURCE_DIR}/libs/arpack-ng/FindArpack.cmake)
    include_directories(${ARPACK_INC_DIR})
    set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES})
    message(STATUS "FOUND PREVIOUSLY INSTALLED ARPACK:   ${ARPACK_LIBRARIES}")
    return()
endif()


message(STATUS "SEARCHING FOR ARPACK IN SYSTEM...")
find_library(ARPACK_LIBRARIES NAMES arpack arpack-ng libarpack libarpack2 libarpack2-dev)
if(ARPACK_LIBRARIES)
    message(STATUS "FOUND ARPACK:   ${ARPACK_LIBRARIES}")
    return()
elseif(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
    message(STATUS "DOWNLOADING ARPACK...")
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
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "ARPACK not found and failed to install: ${res_var}")
    endif()
    set(ARPACK_LIBRARIES "")
    include(${PROJECT_SOURCE_DIR}/libs/arpack-ng/FindArpack.cmake)
    find_package(arpack-ng PATHS "${ARPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
    if(arpack-ng_FOUND)
        get_cmake_property(_variableNames VARIABLES)
        foreach (_variableName ${_variableNames})
            if("${_variableName}" MATCHES "arpack*")
                message(STATUS "${_variableName}=${${_variableName}}")
            endif()
        endforeach()
        include_directories(${ARPACK_INC_DIR})
        set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES})
        message(STATUS "SUCCESSFULLY INSTALLED ARPACK:   ${ARPACK_LIBRARIES}")
        message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/arpack-ng/log_build.txt")
    endif()

else()
    message("Please install LAPACK and BLAS before ARPACK.")
endif()

#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #        if("${_variableName}" MATCHES "Arpack")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #        endif()
#endforeach()