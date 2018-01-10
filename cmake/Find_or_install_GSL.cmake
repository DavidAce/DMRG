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
message("SEARCHING FOR PRE-INSTALLED LIBRARIES: GSL")
find_package(GSL)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    message("FOUND PRE-INSTALLED GSL:   ${GSL_LIBRARIES}")
else()
    message("DOWNLOADING GSL...")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/gsl
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
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
    include(${PROJECT_SOURCE_DIR}/libs/gsl/FindGSL.cmake)
    set(GSL_LIBRARIES "")
    get_libraries(${GSL_LIB_DIR} gsl  GSL_LIBRARIES)
    message("SUCCESSFULLY INSTALLED GSL:   ${GSL_LIBRARIES}")
    message("BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/gsl/log_build.txt")

endif()


