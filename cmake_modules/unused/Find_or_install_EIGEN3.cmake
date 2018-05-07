# The following script will find a native Eigen3 installation on your system.
#
# If these packages are not found, they will be downloaded from a repository and installed
# Either way, this script will define the following variables:
#

# EIGEN3_FOUND - system has eigen lib with correct version
# EIGEN3_INCLUDE_DIR - the eigen include directory
# EIGEN3_VERSION - eigen version
#
# After execution, it is enough to include directory:
#   include_directory(${EIGEN3_INCLUDE_DIR}), which is done here.
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_EIGEN3.cmake)


message("SEARCHING FOR LIBRARY: EIGEN3")
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake")
    include(${PROJECT_SOURCE_DIR}/libs/eigen3/FindEigen3.cmake)
    find_package(Eigen3 3.3 PATHS "${EIGEN3_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE)
    if(EIGEN3_FOUND)
        target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})
        message(STATUS "FOUND PREVIOUSLY INSTALLED EIGEN3: ${EIGEN3_INCLUDE_DIR}")
        return()
    endif()
        message(WARNING "Found 'libs/eigen3/FindEigen3.cmake' in , but failed to configure. You may need a clean build.")
endif()


message(STATUS "SEARCHING FOR EIGEN3 IN SYSTEM...")
find_package(Eigen3 3.3)                                ### Find and define includes for Eigen Library
if(EIGEN3_FOUND)
    message(STATUS "FOUND EIGEN3:   ${EIGEN3_INCLUDE_DIR}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})
    return()
else()

    message(STATUS "DOWNLOADING EIGEN3...")
    set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)


    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/eigen3
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
            -G ${CMAKE_GENERATOR} ../../eigen3
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3)
    execute_process(
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "Eigen3 not found and failed to install: ${res_var}")
    endif()

    # Generate a file that can be detected if the library was installed successfully
    set(EIGEN3_CMAKE_DIR ${INSTALL_DIRECTORY}/eigen3/share/eigen3/cmake)
    file(WRITE ${INSTALL_DIRECTORY}/eigen3/FindEigen3.cmake "set(EIGEN3_CMAKE_DIR     ${EIGEN3_CMAKE_DIR})\n")

    #Include that file
    include(${INSTALL_DIRECTORY}/eigen3/FindEigen3.cmake)
    find_package(Eigen3 3.3 PATHS "${EIGEN3_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
    if(EIGEN3_FOUND)
        target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})
        message(STATUS "SUCCESSFULLY INSTALLED EIGEN3:   ${EIGEN3_INCLUDE_DIR}")
        message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3/log_build.txt")
    endif()
endif()


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#
#
#execute_process(
#        COMMAND ${CMAKE_COMMAND}
#        -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#        -G ${CMAKE_GENERATOR} ../../eigen3
#        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3)
#execute_process(
#        OUTPUT_QUIET
#        OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3/log_build.txt
#        COMMAND ${CMAKE_COMMAND} --build .
#        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/eigen3
#        RESULT_VARIABLE res_var)