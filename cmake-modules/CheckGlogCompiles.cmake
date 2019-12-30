function(check_glog_compiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGETS)
    if(REQUIRED_LIBRARIES_UNPARSED AND REQUIRED_TARGETS)
        message(FATAL_ERROR "Please use EITHER required libs or required targets to honor linking order")
    endif()

    message("CMAKE_VERBOSE_MAKEFILE ${CMAKE_VERBOSE_MAKEFILE}")
    if(CMAKE_VERBOSE_MAKEFILE)
        message(STATUS "Checking if glog compiles")
        message(STATUS "GLOG TEST [${TAG}] REQUIRED_FLAGS                 ${REQUIRED_FLAGS}")
        message(STATUS "GLOG TEST [${TAG}] REQUIRED_DEFINITIONS           ${REQUIRED_DEFINITIONS}")
        message(STATUS "GLOG TEST [${TAG}] REQUIRED_LIBRARIES_UNPARSED    ${REQUIRED_LIBRARIES_UNPARSED}")
        message(STATUS "GLOG TEST [${TAG}] REQUIRED_INCLUDES              ${REQUIRED_INCLUDES}")
        message(STATUS "GLOG TEST [${TAG}] REQUIRED_TARGETS               ${REQUIRED_TARGETS}")
    endif()


    if(NOT BUILD_SHARED_LIBS)
        list(APPEND REQUIRED_LIBRARIES -static)
    endif()

    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES "${expanded_libs}")
            list(APPEND REQUIRED_INCLUDES  "${expanded_incs}")
            list(APPEND REQUIRED_FLAGS     "${expanded_opts}")
        else()
            list(APPEND REQUIRED_LIBRARIES "${elem}")
        endif()
    endforeach()

    foreach(elem ${REQUIRED_TARGETS})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES "${expanded_libs}")
            list(APPEND REQUIRED_INCLUDES  "${expanded_incs}")
            list(APPEND REQUIRED_FLAGS     "${expanded_opts}")
        endif()
    endforeach()
    string(REPLACE ";" " " REQUIRED_FLAGS      "${REQUIRED_FLAGS}") # Needs to be a space-separated list
    list(FILTER REQUIRED_LIBRARIES EXCLUDE REGEX "NOTFOUND") # Make sure every item is valid
    list(FILTER REQUIRED_INCLUDES  EXCLUDE REGEX "NOTFOUND")

    #   Test features
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
    if(CMAKE_VERBOSE_MAKEFILE)
        message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
        message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
    endif()
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
           #include <glog/logging.h>
           int main(int argc, char* argv[]) {
             // Initialize Google's logging library.
             google::InitGoogleLogging(argv[0]);
             LOG(INFO) << \"Found \" << 0 << \" cookies\";
             return 0;
           }
        " GLOG_COMPILES_${TAG})
    if(GLOG_COMPILES_${TAG})
        set(GLOG_COMPILES_${TAG} TRUE PARENT_SCOPE)
    else()
        set(GLOG_COMPILES_${TAG} FALSE PARENT_SCOPE)
        if(CMAKE_VERBOSE_MAKEFILE AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(${ERROR_LOG})
        endif()
    endif()
endfunction()
