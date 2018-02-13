# The following script will find a native GSL installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# GSL_LIBRARIES              - Path to libgsl*.a files
#
# GSL_LIB_DIR                - Path do GSL library directory
# GSL_INC_DIR                - Path to GSL include directory
# GSL_BIN_DIR                - Path to GSL binary directory
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${GSL_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_GSL.cmake)
message("SEARCHING FOR LIBRARY: GSL")
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/gsl/FindGSL.cmake")
    include(${PROJECT_SOURCE_DIR}/libs/gsl/FindGSL.cmake)
    get_libraries(${GSL_LIB_DIR} gsl  GSL_LIBRARIES)
    message(STATUS "FOUND PREVIOUSLY INSTALLED GSL:   ${GSL_LIBRARIES}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
    return()
endif()
message(STATUS "SEARCHING FOR GSL IN SYSTEM...")
find_package(GSL)
if (GSL_FOUND)
    message("FOUND PRE-INSTALLED GSL:   ${GSL_LIBRARIES}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})

else()
    message(STATUS "DOWNLOADING GSL...")
    set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/gsl
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
            -G ${CMAKE_GENERATOR} ../../gsl
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/gsl )
    execute_process(
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/gsl/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/gsl
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "GSL not found and failed to install: ${res_var}")
    endif()
    # Generate a file that can be detected if the library was installed successfully
    set(GSL_LIB_DIR ${INSTALL_DIRECTORY}/gsl/lib)
    set(GSL_INC_DIR ${INSTALL_DIRECTORY}/gsl/include)
    set(GSL_BIN_DIR ${INSTALL_DIRECTORY}/gsl/bin)

    file(WRITE  ${INSTALL_DIRECTORY}/gsl/FindGSL.cmake  "set(GSL_LIB_DIR   ${GSL_LIB_DIR})\n")
    file(APPEND ${INSTALL_DIRECTORY}/gsl/FindGSL.cmake  "set(GSL_INC_DIR   ${GSL_INC_DIR})\n")
    file(APPEND ${INSTALL_DIRECTORY}/gsl/FindGSL.cmake  "set(GSL_BIN_DIR   ${GSL_BIN_DIR})\n")
    file(APPEND ${INSTALL_DIRECTORY}/gsl/FindGSL.cmake  "set(GSL_INCLUDE_DIRS   ${GSL_INC_DIR})\n") # For compatibility with find_package(GSL...)

    #Include that file
    include(${INSTALL_DIRECTORY}/gsl/FindGSL.cmake)
    set(GSL_LIBRARIES "")
    get_libraries(${GSL_LIB_DIR} gsl  GSL_LIBRARIES)
    message(STATUS "SUCCESSFULLY INSTALLED GSL:   ${GSL_LIBRARIES}")
    message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/gsl/log_build.txt")
    target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
endif()


