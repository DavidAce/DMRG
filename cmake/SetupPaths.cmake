cmake_minimum_required(VERSION 3.15)

# Append search paths for find_package and find_library calls
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)


# Transform CMAKE_INSTALL_PREFIX to full path
if(DEFINED CMAKE_INSTALL_PREFIX
    AND NOT IS_ABSOLUTE CMAKE_INSTALL_PREFIX
    AND NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    get_filename_component(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}
            ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR} CACHE FORCE)
    message(STATUS "Setting absolute path CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()

# Setup build and install directories for dependencies
if(NOT DMRG_DEPS_BUILD_DIR)
    set(DMRG_DEPS_BUILD_DIR ${CMAKE_BINARY_DIR}/pkg-build)
endif()

# Install dependencies to the same location as the main project by default
if(NOT DMRG_DEPS_INSTALL_DIR)
    set(DMRG_DEPS_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif()

set(PKG_INSTALL_DIR_DEFAULT ${DMRG_DEPS_INSTALL_DIR} CACHE STRING "" FORCE )
set(PKG_BUILD_DIR_DEFAULT   ${DMRG_DEPS_BUILD_DIR}   CACHE STRING "" FORCE )

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PKG_INSTALL_DIR_DEFAULT} ${CMAKE_INSTALL_PREFIX})
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE INTERNAL "Paths for find_package lookup" FORCE)

if (DMRG_PACKAGE_MANAGER MATCHES "conan")
    # Paths to search for conan installation.
    list(APPEND DMRG_CONAN_HINTS
            ${CONAN_PREFIX}
            $ENV{CONAN_PREFIX}
            ${CONDA_PREFIX}
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3
            $ENV{HOME}/anaconda
            $ENV{HOME}/miniconda3
            $ENV{HOME}/miniconda
            )
    list(APPEND DMRG_CONAN_PATH_SUFFIXES
            bin envs/dmrg/bin
            )
    mark_as_advanced(DMRG_CONAN_HINTS)
    mark_as_advanced(DMRG_CONAN_PATH_SUFFIXES)
endif ()


