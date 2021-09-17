cmake_minimum_required(VERSION 3.18)
project(DMRG++
        VERSION 2.1.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")


# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

# Set options this build
option(DMRG_ENABLE_THREADS "Use C++11 stl threads in Eigen::Tensor" OFF)
option(DMRG_ENABLE_MKL "Enable Intel Math Kernel Library" OFF)
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

# Make an "enum" for valid download methods
set(DMRG_PACKAGE_MANAGERS_VALID find cmake conan)
set(DMRG_PACKAGE_MANAGER find CACHE STRING "Download method for external dependencies")
set_property(CACHE DMRG_PACKAGE_MANAGER PROPERTY STRINGS ${DMRG_PACKAGE_MANAGERS_VALID})
if (NOT DMRG_PACKAGE_MANAGER IN_LIST DMRG_PACKAGE_MANAGERS_VALID)
    message(FATAL_ERROR "DMRG_PACKAGE_MANAGER must be one of ${DMRG_PACKAGE_MANAGERS_VALID}")
endif ()

# Make an "enum" for valid march
set(DMRG_MICROARCH_VALID generic haswell zen zenver1 native)
set(DMRG_MICROARCH native CACHE STRING "CPU micro-architecture")
set_property(CACHE DMRG_MICROARCH PROPERTY STRINGS ${DMRG_MICROARCH_VALID})
if (NOT DMRG_MICROARCH IN_LIST DMRG_MICROARCH_VALID)
    message(FATAL_ERROR "DMRG_MICROARCH must be one of ${DMRG_MICROARCH_VALID}")
endif ()

# Setup paths that find_package should search and
# let cmake find our Find<package>.cmake modules
include(cmake/SetupPaths.cmake)


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake/gitversion.cmake)


############################################################
### Print extra info during CMake configure              ###
############################################################
include(cmake/PrintHostInfo.cmake)
include(cmake/PrintOptionInfo.cmake)


###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
include(cmake/CompilerFlags.cmake)


##################################################
###  Add a target to collect common settings   ###
##################################################
add_library(dmrg-flags INTERFACE)
target_include_directories(dmrg-flags INTERFACE source)
##################################################
###  Add optional RELEASE/DEBUG compile flags  ###
##################################################
# -fsanitize=leak -fsanitize=undefined -fsanitize=alignment  -fsanitize=leak -fsanitize=undefined -fsanitize=alignment
target_compile_options(dmrg-flags INTERFACE -fdiagnostics-color=always)
target_compile_options(dmrg-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wpedantic -Wextra -Wconversion -Wunused>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:-g ${MARCH} ${MTUNE}>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:Clang>>: -fstandalone-debug>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:MINSIZEREL>>:>)
##############################
###  Enable c++17 support  ###
##############################
target_compile_features(dmrg-flags INTERFACE cxx_std_17)

################################
###  Enable build profiling  ###
################################
if (DMRG_PROFILE_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(dmrg-flags INTERFACE -ftime-trace)
endif ()

#######################################
# Settings for sanitizers           ###
#######################################
if (DMRG_ENABLE_ASAN)
    target_compile_options(dmrg-flags INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_libraries(dmrg-flags INTERFACE -fsanitize=address)
    #    if(NOT BUILD_SHARED_LIBS)
    #        target_link_libraries(dmrg-flags INTERFACE -static-libasan)
    #    endif()
endif ()
if (DMRG_ENABLE_USAN)
    target_compile_options(dmrg-flags INTERFACE -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_libraries(dmrg-flags INTERFACE -fsanitize=undefined)
endif ()

##################################################################
###  Link system libs statically                               ###
##################################################################
if (NOT BUILD_SHARED_LIBS)
    target_link_options(dmrg-flags INTERFACE -static-libgcc -static-libstdc++)
endif ()


##################################################################
### Speed up compilation with precompiled headers              ###
##################################################################
if (DMRG_ENABLE_PCH)
    add_library(dmrg-pch OBJECT cmake/empty.cpp)
    target_link_libraries(dmrg-pch PRIVATE dmrg-deps dmrg-flags)
    target_precompile_headers(dmrg-pch PUBLIC
            <string> <string_view> <vector> <array>
            <optional> <complex> <memory> <chrono> <algorithm>
            <fmt/compile.h>
            <fmt/core.h>
            <fmt/format.h>
            <fmt/ostream.h>
            <fmt/ranges.h>
            <spdlog/spdlog.h>
            <spdlog/sinks/stdout_color_sinks.h>
            <Eigen/Core>
            <unsupported/Eigen/CXX11/Tensor>
            <h5pp/h5pp.h>
            )
endif ()

#######################################
# Settings for ccache               ###
#######################################
if (DMRG_ENABLE_CCACHE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0.0")
        set(COMPILER_OK FALSE)
    else ()
        set(COMPILER_OK TRUE)
    endif ()
    mark_as_advanced(COMPILER_OK)
    if (COMPILER_OK)
        find_program(CMAKE_CXX_COMPILER_LAUNCHER ccache REQUIRED)
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(dmrg-pch PUBLIC -Xclang -fno-pch-timestamp -fpch-preprocess)
        endif ()
        message(STATUS "Using ccache ${CMAKE_CXX_COMPILER_LAUNCHER}")
        if (DMRG_ENABLE_PCH)
            message(STATUS "Detected ccache + pch: Remember to set --> sloppiness = include_file_mtime,pch_defines,time_macros <-- in your ccache.conf")
        endif ()
    endif ()
endif ()


###


#######################################
### Compiler-dependent linker flags ###
#######################################
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(dmrg-flags INTERFACE -stdlib=libstdc++)
endif ()


###########################################
### Fetch or find all the dependencies  ###
###########################################
include(cmake/SetupDependencies.cmake)


################################
###  Add all source files    ###
################################
add_executable(DMRG++)
set_target_properties(DMRG++ PROPERTIES LINK_WHAT_YOU_USE TRUE)


##################################################################
### Check if we can enable Link Time Optimization              ###
##################################################################
if(DMRG_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(DMRG++ PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
        if(CMAKE_VERSION VERSION_GREATER 3.19)
            set_target_properties(dmrg-flags PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
        endif()
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
        set(DMRG_ENABLE_LTO OFF CACHE)
    endif ()
endif ()

add_subdirectory(source)

#######################################################
###  Check  #include<optional> or                   ###
###  or     #include<experimental/optional>         ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake/CheckCXXOptional.cmake)
CheckCXXOptional()


#######################################


##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(DMRG++ PRIVATE dmrg-main)  # Link everything to main project
#target_sources(DMRG++ PRIVATE $<TARGET_OBJECTS:dmrg-main>)


# Uninstall target
if(NOT TARGET uninstall)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DMRGUninstall.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake
            IMMEDIATE @ONLY)

    add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake)
endif()



# Print summary of CMake configuration
include(cmake/PrintTargetInfo.cmake)
include(cmake/getExpandedTarget.cmake)
message(DEBUG "| DMRG TARGET SUMMARY")
message(DEBUG "|-----------------------------")
print_compiler_info("| ")
expand_target_all_targets(DMRG++ DMRG_TARGET_EXPANDED)
foreach (tgt ${DMRG_TARGET_EXPANDED})
    print_target_info(${tgt} "| ")
endforeach ()


###################
### INFORMATION ###
###################

message("")
message("=============================================================================")
message("")
if (OPENBLAS_DYNAMIC_ARCH AND NOT TARGET mkl::mkl)
    # Using OpenBLAS
    message("   OpenBLAS DYNAMIC_ARCH: Remember to set architecture before launch:")
    message("          export OPENBLAS_CORETYPE=<microarch>\n")
endif ()
message("   To clean cmake files, run './build.sh -c'")
message("   To see more options,  run './build.sh -h'")
message("")
message("=============================================================================")
message("")


#############################
## UNIT TEST
#############################
# Simple testing with ctest
if (DMRG_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()


#############################
## DOCS
#############################
if (DMRG_ENABLE_DOCS)
    add_subdirectory(docs)
endif ()






