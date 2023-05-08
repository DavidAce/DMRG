cmake_minimum_required(VERSION 3.15)
project(external-mpreal)
include(ExternalProject)
ExternalProject_Add(external_mpreal
        URL         https://github.com/advanpix/mpreal/archive/refs/tags/mpfrc++-3.6.9.tar.gz
#        URL         https://github.com/advanpix/mpreal/archive/refs/tags/mpfrc++-3.6.8.tar.gz
        URL_MD5     80a9c8665b000eab85b2342fa0951b2b # mpfrc++-3.6.9.tar.gz
#        URL_MD5     737311686cedb7e3da36e1c2b5edc5df # mpfrc++-3.6.8.tar.gz
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS TRUE
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
        ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/mpreal.h <INSTALL_DIR>/include/mpreal.h
        )