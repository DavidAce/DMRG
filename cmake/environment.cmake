cmake_minimum_required(VERSION 3.24)
function(scan_environment)
    #######################
    ### Get git version ###
    #######################
    execute_process(
            COMMAND git rev-parse
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE result_var
            OUTPUT_VARIABLE GIT_REV_PARSE
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(result_var STREQUAL "0")
        # Get the current working branch
        execute_process(
                COMMAND git rev-parse --abbrev-ref HEAD
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_BRANCH
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Get the latest abbreviated commit hash of the working branch
        execute_process(
                COMMAND git rev-parse --verify HEAD
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Get the revision count
        execute_process(
                COMMAND git rev-list HEAD --count
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_REVISION
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
    cmake_host_system_information(RESULT _proc_type QUERY PROCESSOR_DESCRIPTION)
    cmake_host_system_information(RESULT _os_name QUERY OS_NAME)
    cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
    cmake_host_system_information(RESULT _os_version QUERY OS_VERSION)
    cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)

    configure_file(
            ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/environment.h.in
            ${CMAKE_BINARY_DIR}/environment/include/env/environment.h
    )

    include_directories(${CMAKE_BINARY_DIR}/environment/include)
endfunction()

scan_environment()

