cmake_minimum_required(VERSION 3.15)

# Append search paths for find_package and find_library calls
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)


# Transform CMAKE_INSTALL_PREFIX to full path
if(DEFINED CMAKE_INSTALL_PREFIX AND NOT IS_ABSOLUTE CMAKE_INSTALL_PREFIX)
    get_filename_component(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}
            ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR} CACHE FORCE)
    message(STATUS "Setting absolute path CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()

# Setup build and install directories for dependencies
if(NOT DMRG_DEPS_BUILD_DIR)
    set(DMRG_DEPS_BUILD_DIR ${CMAKE_BINARY_DIR}/dmrg-deps-build)
endif()
if(NOT DMRG_DEPS_INSTALL_DIR)
    set(DMRG_DEPS_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}) # Install to the same location as h5pp by default
endif()

if(DMRG_PREFIX_ADD_PKGNAME)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/dmrg CACHE STRING
            "The option DMRG_PREFIX_ADD_PKGNAME=ON sets the install directory: <CMAKE_INSTALL_PREFIX>/dmrg"
            FORCE)
endif()

# Add search directories and flags for the CMake find_* tools
if(DMRG_PACKAGE_MANAGER STREQUAL "find")
    set(REQUIRED REQUIRED)
endif()
if(DMRG_PACKAGE_MANAGER STREQUAL "cmake")
    # We set variables here that allows us to find packages exclusively with CMAKE_PREFIX_PATH
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
    # Flags that can be used directly on find_package
    # Enumerated according to the cmake manual for find_package
    set(N5 NO_SYSTEM_ENVIRONMENT_PATH) #5
    set(N6 NO_CMAKE_PACKAGE_REGISTRY) #6
    set(N7 NO_CMAKE_SYSTEM_PATH) #7
    set(N8 NO_CMAKE_SYSTEM_PACKAGE_REGISTRY) #8

endif()
if(DMRG_PACKAGE_MANAGER MATCHES "cmake")
    list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX})
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE STRING "" FORCE)
endif()



# Make sure find_library prefers static/shared library depending on BUILD_SHARED_LIBS
# This is important when finding dependencies such as zlib which provides both shared and static libraries.
# Note that we do not force this cache variable, so users can override it
#if(BUILD_SHARED_LIBS)
#    # This is order is the default
#    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX};${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE STRING "Prefer finding shared libraries")
#else()
#    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX};${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE STRING "Prefer finding static libraries")
#endif()

if (CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES "64")
    set(FIND_LIBRARY_USE_LIB64_PATHS ON)
elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(FIND_LIBRARY_USE_LIB32_PATHS ON)
endif ()


if(WIN32)
    # On Windows it is standard practice to collect binaries into one directory.
    # This way we avoid errors from .dll's not being found at runtime.
    # These directories will contain h5pp tests, examples and possibly dependencies
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" CACHE PATH "Collect .exe and .dll")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" CACHE PATH "Collect .lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" CACHE PATH "Collect .lib")
endif()

if (DMRG_PACKAGE_MANAGER MATCHES "conan")
    # Paths to search for conan installation.
    list(APPEND DMRG_CONAN_CANDIDATE_PATHS
            ${CONAN_PREFIX}
            $ENV{CONAN_PREFIX}
            ${CONDA_PREFIX}
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3/envs/dmrg
            $ENV{HOME}/anaconda/envs/dmrg
            $ENV{HOME}/miniconda3/envs/dmrg
            $ENV{HOME}/miniconda/envs/dmrg
            $ENV{HOME}/.conda/envs/dmrg
            $ENV{HOME}/anaconda3
            $ENV{HOME}/anaconda
            $ENV{HOME}/miniconda3
            $ENV{HOME}/miniconda
            $ENV{HOME}/.conda
            )
endif ()


mark_as_advanced(DMRG_CONAN_CANDIDATE_PATHS)