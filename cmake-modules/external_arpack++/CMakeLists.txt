cmake_minimum_required(VERSION 3.12)
project(external-arpack++)

include(ExternalProject)
ExternalProject_Add(external_ARPACK++
        URL         https://github.com/m-reuter/arpackpp/archive/2.3.0.tar.gz
        URL_MD5     1b09e35b6c44e118003922643b99978a
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND ""
        INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
        )