cmake_minimum_required(VERSION 3.24)
project(DMRG++
        VERSION 3.1.0
        DESCRIPTION "MPS algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")

include(cmake/Deprecation.cmake) # Warn if user passes deprecated settings

# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/modules)

# Set options for this build
option(DMRG_USE_QUADMATH "Use __float128 instead of long double for fLBIT time evolution" FALSE)
option(DMRG_ENABLE_TBLIS "Use TBLIS for real-valued tensor contractions instead of Eigen3" FALSE)
option(DMRG_ENABLE_TESTS "Enable unit testing with ctest" FALSE)
option(DMRG_BUILD_EXAMPLES "Build examples" FALSE)
option(DMRG_BUILD_TOOLS "Build tools" FALSE)
option(DMRG_ENABLE_DOCS "Build documentation" FALSE)
option(DMRG_CMAKE_DEBUG "Extra information during CMake configuration" FALSE)
option(CMAKE_POSITION_INDEPENDENT_CODE "Use -fPIC when compiling shared libraries" TRUE)
option(COMPILER_PROFILE_BUILD "Enable -ftime-trace (inspect with ClangBuildAnalyzer)" FALSE)
option(COMPILER_ENABLE_ASAN "Enable runtime address sanitizer -fsanitize=address" FALSE)
option(COMPILER_ENABLE_USAN "Enable undefined behavior sanitizer -fsanitize=undefined" FALSE)
option(COMPILER_ENABLE_LTO "Enable link time optimization" FALSE)
option(COMPILER_ENABLE_PCH "Enable precompiled headers to speed up compilation" FALSE)
option(COMPILER_ENABLE_CCACHE "Enable ccache to speed up compilation" FALSE)
option(COMPILER_ENABLE_COVERAGE "Enable test coverage" FALSE)
option(EIGEN_USE_THREADS "Use STL threads to parallelize Eigen::Tensor" TRUE)
option(THREADS_PREFER_PTHREAD_FLAG "Prefer -pthread flag over -lpthread or similar" TRUE)
option(CMAKE_INTERPROCEDURAL_OPTIMIZATION "Enable link-time-optimization" FALSE)
option(DMRG_BENCH_CONTRACTION "For internal development"  FALSE)
option(DMRG_SAVE_CONTRACTION "For internal development"  FALSE)


# Make an "enum" for valid package managers
set(DMRG_PACKAGE_MANAGERS_VALID find cmake conan)
set(DMRG_PACKAGE_MANAGER find CACHE STRING "Package manager for external dependencies")
set_property(CACHE DMRG_PACKAGE_MANAGER PROPERTY STRINGS ${DMRG_PACKAGE_MANAGERS_VALID})
if(NOT DMRG_PACKAGE_MANAGER IN_LIST DMRG_PACKAGE_MANAGERS_VALID)
    message(FATAL_ERROR "DMRG_PACKAGE_MANAGER must be one of ${DMRG_PACKAGE_MANAGERS_VALID}")
endif()



###  Check <optional>
include(cmake/CheckCXXOptional.cmake)
checkcxxoptional()

### Print extra info during CMake configure
include(cmake/PrintBuildInfo.cmake)

###  Add targets to collect common settings
add_library(dmrg-flags INTERFACE)
add_library(dmrg-deps INTERFACE)
target_include_directories(dmrg-flags INTERFACE source)

###  Apply compiler flags
include(cmake/CompilerFlags.cmake)

################################################################
### Get environment, host, build, exec and git details       ###
### Generates a header env/environment.h                     ###
### Include it using #include <env/environment.h>            ###
################################################################
include(cmake/environment.cmake)

### Create the main executable
add_executable(DMRG++)
set_target_properties(DMRG++ PROPERTIES LINK_WHAT_YOU_USE TRUE)

### Link precompiled headers. Requires DMRG_ENABLE_PCH=ON
target_link_precompiled_headers(DMRG++)

# Setup paths that find_package should search and
# let cmake find our Find<package>.cmake modules
include(cmake/SetupPaths.cmake)

### Find or install all the dependencies
include(cmake/SetupDependencies.cmake)

###  Add all source files
add_subdirectory(source)

### Link all the things!
target_link_libraries(DMRG++ PRIVATE dmrg-main)

### Enable static linking of libgcc and libstdc++. Requires BUILD_SHARED_LIBS=OFF
target_enable_static_libgcc(DMRG++)

# Uninstall target
if(NOT TARGET uninstall)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DMRGUninstall.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake
            IMMEDIATE @ONLY)

    add_custom_target(uninstall
                      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake)
endif()

if(DMRG_CMAKE_DEBUG)
    # Print summary of CMake configuration
    include(cmake/PrintTargetInfo.cmake)
    print_and_write_project_summary(DMRG++)
endif()

if(DMRG_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests) # Unit testing with ctest and catch2
endif()
add_subdirectory(bench) # Benchmarks with ctest and nanobench
add_subdirectory(tools) # Build tools for minor DMRG related tasks
add_subdirectory(docs)  # Build documentation with doxygen

#############################
## NINJA TRACING
#############################
if(CMAKE_GENERATOR MATCHES Ninja)
    add_custom_command(
            TARGET DMRG++
            COMMAND ninja -C ${CMAKE_BINARY_DIR} -t recompact
            COMMAND python ${CMAKE_CURRENT_LIST_DIR}/ninjatracing --showall ${CMAKE_BINARY_DIR}/.ninja_log > ${CMAKE_BINARY_DIR}/.ninja_trace.json
            DEPENDS ${CMAKE_BINARY_DIR}/.ninja_log
            BYPRODUCTS ${CMAKE_BINARY_DIR}/.ninja_trace.json
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            COMMENT "Generating ninja build statistics"
    )
endif()


