
#######################
### Get git version ###
#######################

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

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
cmake_host_system_information(RESULT _proc_type QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT _os_name QUERY OS_NAME)
cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
cmake_host_system_information(RESULT _os_version QUERY OS_VERSION)
cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)


configure_file(
${CMAKE_SOURCE_DIR}/cmake/environment.h.in
${CMAKE_BINARY_DIR}/environment/include/env/environment.h
)

include_directories(${CMAKE_BINARY_DIR}/environment/include)

