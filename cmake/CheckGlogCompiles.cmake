function(check_glog_compiles TARGETS LIBS INCS OPTS DEFS)
    if(NOT BUILD_SHARED_LIBS)
        list(APPEND CMAKE_REQUIRED_LIBRARIES -static)
    endif()
    list(APPEND CMAKE_REQUIRED_LIBRARIES     ${LIBS} ${TARGETS})
    list(APPEND CMAKE_REQUIRED_INCLUDES      ${INCS})
    list(APPEND CMAKE_REQUIRED_FLAGS         ${OPTS})
    list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${DEFS})

    list(TRANSFORM "CMAKE_REQUIRED_DEFINITIONS" PREPEND "-D")  # Definitions should start with "-D"
    string(REPLACE ";" " "  CMAKE_REQUIRED_FLAGS  "${CMAKE_REQUIRED_FLAGS}")        # Needs to be a space-separated list

    if(DMRG_PRINT_CHECKS)
        message(STATUS "GLOG COMPILE TEST CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "GLOG COMPILE TEST CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "GLOG COMPILE TEST CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "GLOG COMPILE TEST CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    endif()
    #   Test features
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
           #include <glog/logging.h>
           int main(int argc, char* argv[]) {
             // Initialize Google's logging library.
             google::InitGoogleLogging(argv[0]);
             LOG(INFO) << \"Found \" << 0 << \" cookies\";
             return 0;
           }
        " GLOG_COMPILES)
    if(NOT GLOG_COMPILES)
        unset(GLOG_COMPILES CACHE)
        unset(GLOG_COMPILES PARENT_SCOPE)
        if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
        endif()
    endif()
endfunction()
