cmake_minimum_required(VERSION 3.15)

set(PROJECT_UNAME DMRG)
set(PROJECT_LNAME dmrg)


message(STATUS "C compiler ${CMAKE_C_COMPILER}")
message(STATUS "FC compiler ${CMAKE_Fortran_COMPILER}")
message(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}")

#####################################################
### Set the  microarchitecture for OpenBLAS       ###
#####################################################
# Make an "enum" for valid march
set(${PROJECT_UNAME}_MICROARCH_VALID generic haswell zen zenver1 native)
set(${PROJECT_UNAME}_MICROARCH native CACHE STRING "CPU micro-architecture")
set_property(CACHE ${PROJECT_UNAME}_MICROARCH PROPERTY STRINGS ${${PROJECT_UNAME}_MICROARCH_VALID})
if (NOT ${PROJECT_UNAME}_MICROARCH IN_LIST ${PROJECT_UNAME}_MICROARCH_VALID)
    message(FATAL_ERROR "${PROJECT_UNAME}_MICROARCH must be one of ${${PROJECT_UNAME}_MICROARCH_VALID}")
endif ()


cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
set(OPENBLAS_TARGET HASWELL)
set(OPENBLAS_DYNAMIC_ARCH ON)
if($ENV{CI} OR $ENV{GITHUB_ACTIONS} OR ${PROJECT_UNAME}_MICROARCH MATCHES "generic")
    set(MARCH -march=x86-64)
    set(MTUNE -mtune=generic)
    set(OPENBLAS_TARGET GENERIC)
    set(OPENBLAS_DYNAMIC_ARCH OFF)
elseif(DEFINED ${PROJECT_UNAME}_MICROARCH)
    set(MARCH -march=${${PROJECT_UNAME}_MICROARCH})
    set(MTUNE -mtune=${${PROJECT_UNAME}_MICROARCH})
else()
    set(MARCH -march=haswell)
    set(MTUNE -mtune=native)
endif()
message(DEBUG "Using ${MARCH} ${MTUNE}")



######################################################################################
###                   Apply RELEASE/DEBUG compile flags                            ###
######################################################################################
# I have benchmarked the compiler flags below
#        -fstack-protector
#        -D_FORTIFY_SOURCE=2
#        -fno-omit-frame-pointer
#        -fno-strict-aliasing
# and found that there is NO significant performance difference in tensor contractions.
# These were the results from Tensorbench, best of 2 runs,
# variability between runs ~0.3 seconds.
# (all with -O3 -mfma -DNDEBUG -march=native -mtune=native):
# 19.3789 s: -fstack-protector
# 19.2768 s: -fstack-protector -fno-omit-frame-pointer -fno-strict-aliasing
# 19.6200 s: -fstack-protector -fno-omit-frame-pointer -fno-strict-aliasing -D_FORTIFY_SOURCE=2
# 19.3540 s: (none)
#
# In particular, -fno-strict-aliasing fixes a bug in tensor shuffle using -O3 -DNDEBUG
# in gcc 10.0.1. The bug causes the dimensions to blow up when assigning the shuffling
# operation to a new tensor, such as in:
#       Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
# Since the dimensions become kind of random (and not near the limits of int/long),
# I believe gcc-10 is simply reordering something it shouldn't, thus assigning
# uninitialized values on these dimensions.
# The bug was hard to track down, and could be fixed by various other flags that
# hurt performance more, like -mno-avx, removing -DNDEBUG, or lowering -O3 to -O2.
######################################################################################


set(CMAKE_EXPORT_COMPILE_COMMANDS ON) ### Write compile commands to file
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -fno-strict-aliasing -fdiagnostics-color=always -Wall -Wpedantic -Wextra -Wconversion -Wunused")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MARCH} ${MTUNE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fstack-protector -D_FORTIFY_SOURCE=2") #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fstack-protector -D_FORTIFY_SOURCE=2")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_EXE_LINKER_FLAGS MATCHES "fuse-ld=gold")
    set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=gold -Wl,--disable-new-dtags")
endif()

# Set these variables so that the same flags are used for building dependencies
set(CMAKE_CXX_FLAGS_INIT                 "${CMAKE_CXX_FLAGS_INIT} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE_INIT         "${CMAKE_CXX_FLAGS_RELEASE_INIT} ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG_INIT           "${CMAKE_CXX_FLAGS_DEBUG_INIT} ${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_INIT          "${CMAKE_EXE_LINKER_FLAGS_INIT} ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)


###############################
# Settings for shared builds
###############################

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when installing)
# Note: Since DMRG++ is often run from the build folder we want to keep the build-folder RPATH in the executable.
#       Therefore it makes sense to keep this setting "FALSE" here but "TRUE" for dependencies that are
#       installed with in "cmake" mode with externalproject_add
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


###  Add optional RELEASE/DEBUG compile to flags
if(NOT TARGET flags)
    add_library(flags INTERFACE)
endif()
target_compile_options(flags INTERFACE -g -fno-strict-aliasing -fdiagnostics-color=always -Wall -Wpedantic -Wextra -Wconversion -Wunused)
target_compile_options(flags INTERFACE $<$<CONFIG:RELEASE>:${MARCH} ${MTUNE}>)
target_compile_options(flags INTERFACE $<$<CONFIG:DEBUG>: -fno-omit-frame-pointer -fstack-protector -D_FORTIFY_SOURCE=2>)
target_compile_options(flags INTERFACE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:Clang>>: -fstandalone-debug>)
target_compile_options(flags INTERFACE $<$<CONFIG:RELWITHDEBINFO>:-fno-omit-frame-pointer -fstack-protector -D_FORTIFY_SOURCE=2>)
target_compile_options(flags INTERFACE $<$<CONFIG:MINSIZEREL>:>)

###  Enable c++17 support
target_compile_features(flags INTERFACE cxx_std_17)

###  Enable build profiling with ClangBuildAnalyzer
if (${PROJECT_UNAME}_PROFILE_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(flags INTERFACE -ftime-trace)
endif ()

# Settings for sanitizers
if (${PROJECT_UNAME}_ENABLE_ASAN)
    target_compile_options(flags INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_libraries(flags INTERFACE -fsanitize=address)
    #    if(NOT BUILD_SHARED_LIBS)
    #        target_link_libraries(flags INTERFACE -static-libasan)
    #    endif()
endif ()
if (${PROJECT_UNAME}_ENABLE_USAN)
    target_compile_options(flags INTERFACE -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_libraries(flags INTERFACE -fsanitize=undefined)
endif ()


###  Link system libs statically
if (NOT BUILD_SHARED_LIBS)
    target_link_options(flags INTERFACE -static-libgcc -static-libstdc++)
endif ()


### Speed up compilation with precompiled headers
function(target_link_precompiled_headers tgt)
    if(${PROJECT_UNAME}_ENABLE_PCH)
        set(PCH_HEADERS
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

        get_target_property(type ${tgt} TYPE)
        if(type MATCHES "EXECUTABLE")
            if(NOT TARGET pch-exe)
                add_executable(pch-exe ${PROJECT_SOURCE_DIR}/cmake/pch.cpp)
                target_link_libraries(pch-exe PUBLIC deps flags)
                target_precompile_headers(pch-exe PUBLIC ${PCH_HEADERS})
            endif()
            target_precompile_headers(${tgt} REUSE_FROM pch-exe)
        else()
            if(NOT TARGET pch-obj)
                add_library(pch-obj OBJECT ${PROJECT_SOURCE_DIR}/cmake/pch.cpp)
                target_link_libraries(pch-obj PUBLIC deps flags)
                target_precompile_headers(pch-obj PUBLIC ${PCH_HEADERS})
            endif()
            target_precompile_headers(${tgt} REUSE_FROM pch-obj)
        endif()
        message(TRACE "Enabled precompiled headers for target ${tgt}")
    endif()
endfunction()


# Settings for ccache
if (${PROJECT_UNAME}_ENABLE_CCACHE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0.0")
        set(COMPILER_OK FALSE)
    else ()
        set(COMPILER_OK TRUE)
    endif ()
    mark_as_advanced(COMPILER_OK)
    if (COMPILER_OK)
        find_program(CMAKE_CXX_COMPILER_LAUNCHER ccache REQUIRED)
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(flags INTERFACE -Xclang -fno-pch-timestamp -fpch-preprocess)
        endif ()
        message(STATUS "Using ccache ${CMAKE_CXX_COMPILER_LAUNCHER}")
        if (${PROJECT_UNAME}_ENABLE_PCH)
            message(STATUS "Detected ccache + pch: Remember to set --> sloppiness = include_file_mtime,pch_defines,time_macros <-- in your ccache.conf")
        endif ()
    endif ()
endif ()


# Compiler-dependent linker flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(flags INTERFACE -stdlib=libstdc++)
endif ()

