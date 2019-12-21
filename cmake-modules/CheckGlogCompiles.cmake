function(check_glog_compiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGET)
    #    message(STATUS "Checking if arpack headers are working")
    set(REQUIRED_LIBRARIES)
    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
        else()
            list(APPEND REQUIRED_LIBRARIES ${elem})
        endif()
    endforeach()

    foreach(elem ${REQUIRED_TARGET})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
            list(APPEND REQUIRED_INCLUDES  ${expanded_incs})
            list(APPEND REQUIRED_FLAGS     ${expanded_opts})
        endif()
    endforeach()
    string(REPLACE ";" " " REQUIRED_FLAGS      "${REQUIRED_FLAGS}") # Needs to be a space-separated list
    #   Test features
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
    message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
    message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
    message(STATUS "GLOG TEST [${TAG}] CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")

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
    endif()
endfunction()
