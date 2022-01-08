cmake_minimum_required(VERSION 3.18)
project(DMRG++
        VERSION 2.1.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")

# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

# Set options this build
option(DMRG_ENABLE_MKL "Enable Intel Math Kernel Library" OFF)
option(DMRG_ENABLE_RSVD "Enable Randomized SVD module" OFF)
option(DMRG_ENABLE_TESTS "Enable unit testing with ctest" OFF)
option(DMRG_ENABLE_ASAN "Enable runtime address sanitizer -fsanitize=address" OFF)
option(DMRG_ENABLE_USAN "Enable undefined behavior sanitizer -fsanitize=undefined" OFF)
option(DMRG_ENABLE_LTO "Enable link time optimization" OFF)
option(DMRG_ENABLE_PCH "Enable precompiled headers to speed up compilation" ON)
option(DMRG_ENABLE_CCACHE "Enable ccache to speed up compilation" ON)
option(DMRG_ENABLE_COVERAGE "Enable test coverage" OFF)
option(DMRG_PROFILE_BUILD "Enable build-time profiling with -ftime-trace (study results with ClangBuildAnalyzer)" OFF)
option(DMRG_BUILD_EXAMPLES "Build examples" OFF)
option(DMRG_DEPS_IN_SUBDIR "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>" ON)
option(DMRG_ENABLE_DOCS "Build documentation" OFF)
option(CMAKE_POSITION_INDEPENDENT_CODE  "Use -fPIC when compiling shared libraries" ON)

# Make an "enum" for valid package managers
set(DMRG_PACKAGE_MANAGERS_VALID find cmake conan)
set(DMRG_PACKAGE_MANAGER find CACHE STRING "Package manager for external dependencies")
set_property(CACHE DMRG_PACKAGE_MANAGER PROPERTY STRINGS ${DMRG_PACKAGE_MANAGERS_VALID})
if (NOT DMRG_PACKAGE_MANAGER IN_LIST DMRG_PACKAGE_MANAGERS_VALID)
    message(FATAL_ERROR "DMRG_PACKAGE_MANAGER must be one of ${DMRG_PACKAGE_MANAGERS_VALID}")
endif ()

###  Check <optional>
include(cmake/CheckCXXOptional.cmake)
CheckCXXOptional()


################################################################
### Get environment, host and git details                    ###
### Generates a header env/environment.h                     ###
### Include it using #include <env/environment.h>            ###
################################################################
include(cmake/environment.cmake)


### Print extra info during CMake configure
include(cmake/PrintBuildInfo.cmake)


###  Add targets to collect common settings
add_library(deps INTERFACE)
add_library(flags INTERFACE)
target_include_directories(flags INTERFACE source)

###  Apply compiler flags
include(cmake/CompilerFlags.cmake)


# Setup paths that find_package should search and
# let cmake find our Find<package>.cmake modules
include(cmake/SetupPaths.cmake)


### Find or install all the dependencies
include(cmake/SetupDependencies.cmake)



###  Add all source files
add_executable(DMRG++)
set_target_properties(DMRG++ PROPERTIES LINK_WHAT_YOU_USE TRUE)
add_subdirectory(source)

### Link all the things!
target_link_libraries(DMRG++ PRIVATE dmrg-main)

### Link precompiled headers. Requires DMRG_ENABLE_PCH=ON
target_link_precompiled_headers(DMRG++)

### Enable link time optimization. Requires DMRG_ENABLE_LTO=ON
target_enable_lto(DMRG++)


# Uninstall target
if (NOT TARGET uninstall)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DMRGUninstall.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake
            IMMEDIATE @ONLY)

    add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake)
endif()



# Print summary of CMake configuration
include(cmake/PrintTargetInfo.cmake)
print_project_summary(DMRG++)


###################
### INFORMATION ###
###################

message(STATUS "")
message(STATUS "=============================================================================")
message(STATUS "")
if (OPENBLAS_DYNAMIC_ARCH AND NOT TARGET mkl::mkl)
    # Using OpenBLAS
    message(STATUS "   OpenBLAS DYNAMIC_ARCH: Remember to set architecture before launch:")
    message(STATUS "          export OPENBLAS_CORETYPE=<microarch>\n")
endif ()
message(STATUS "   To clean cmake files, run './build.sh -c'")
message(STATUS "   To see more options,  run './build.sh -h'")
message(STATUS "")
message(STATUS "=============================================================================")
message(STATUS "")


#############################
## UNIT TEST
#############################
# Simple testing with ctest
if (DMRG_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

#############################
## BENCHMARKS
#############################
# Simple testing with ctest
if (DMRG_ENABLE_BENCH)
    enable_testing()
    add_subdirectory(bench)
endif ()


#############################
## DOCS
#############################
if (DMRG_ENABLE_DOCS)
    add_subdirectory(docs)
endif ()






