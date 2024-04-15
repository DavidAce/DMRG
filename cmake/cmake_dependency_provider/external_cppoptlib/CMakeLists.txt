cmake_minimum_required(VERSION 3.15)
project(external-cppoptlib)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(ExternalProject)
ExternalProject_Add(external_cppoptlib
        GIT_REPOSITORY  https://github.com/PatWie/CppNumericalSolvers.git
        GIT_TAG 9e21e01736237a1e10e57bffca1b00fe6655d9f6
        GIT_SHALLOW TRUE
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}

        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        #Give Ninja access to the terminal.
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_UPDATE TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
        USES_TERMINAL_TEST TRUE

        UPDATE_COMMAND ""
        TEST_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND ""
        INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include

        )