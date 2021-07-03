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
option(DMRG_BUILD_EXAMPLES "Build examples" OFF)
option(DMRG_PRINT_INFO "Print information during cmake configure" OFF)
option(DMRG_PRINT_CHECKS "Print more information during cmake configure" OFF)
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


################################
###  Add all source files    ###
################################
include(source_list.txt)
add_executable(DMRG++)
set_target_properties(DMRG++ PROPERTIES LINK_WHAT_YOU_USE TRUE)
add_library(dmrg-main OBJECT ${DMRG_SOURCES})

add_library(dmrg-arp OBJECT
        # This one is special, since it absolutely can't have BLAS/LAPACK headers during compile for some reason
        source/math/eig/arpack_solver/arpack_solver.cpp
        )

add_library(dmrg-eig OBJECT
        source/math/eig/settings.cpp
        source/math/eig/solver.cpp
        source/math/eig/solver_dgeev.cpp
        source/math/eig/solver_dsyevd.cpp
        source/math/eig/solver_zgeev.cpp
        source/math/eig/solver_zheevd.cpp
        source/math/eig/solver_primme.cpp
        source/math/eig/solution.cpp
        source/math/eig/matvec/matvec_dense.cpp
        source/math/eig/matvec/matvec_sparse.cpp
        source/math/eig/matvec/matvec_mpo.cpp
        )

add_library(dmrg-opt OBJECT
        source/tools/finite/opt_mps.cpp
        source/tools/finite/opt-internal/ceres_base.cpp
        source/tools/finite/opt-internal/ceres_direct.cpp
        source/tools/finite/opt-internal/ceres_direct_functor.cpp
        source/tools/finite/opt-internal/ceres_subspace.cpp
        source/tools/finite/opt-internal/ceres_subspace_diag_full.cpp
        source/tools/finite/opt-internal/ceres_subspace_diag_part.cpp
        source/tools/finite/opt-internal/ceres_subspace_filter.cpp
        source/tools/finite/opt-internal/ceres_subspace_functor.cpp
        source/tools/finite/opt-internal/ceres_subspace_hamiltonian.cpp
        source/tools/finite/opt-internal/ground_state.cpp
        source/tools/finite/opt-internal/krylov_energy.cpp
        source/tools/finite/opt-internal/krylov_variance.cpp
        source/tools/finite/opt-internal/krylov_solution.cpp
        source/tools/finite/opt-internal/report.cpp
        source/tools/finite/opt.cpp
        )

add_library(dmrg-dbg OBJECT
        source/general/stack_trace.cpp
        )

add_library(dmrg-primme INTERFACE)

target_include_directories(dmrg-main PUBLIC source)
target_include_directories(dmrg-eig PUBLIC source)
target_include_directories(dmrg-arp PUBLIC source)
target_include_directories(dmrg-opt PUBLIC source)
target_include_directories(dmrg-dbg PUBLIC source)

add_library(dmrg-flags INTERFACE)
target_include_directories(dmrg-flags SYSTEM INTERFACE source)

add_library(dmrg-deps INTERFACE)

#######################################################
###  Check  #include<optional> or                   ###
###  or     #include<experimental/optional>         ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake/CheckCXXOptional.cmake)
CheckCXXOptional()



##################################################
###  Add optional RELEASE/DEBUG compile flags  ###
##################################################
# -fsanitize=leak -fsanitize=undefined -fsanitize=alignment  -fsanitize=leak -fsanitize=undefined -fsanitize=alignment
target_compile_options(dmrg-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:-g ${MARCH} ${MTUNE}>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>: -Wpedantic -Wextra -Wconversion>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:MINSIZEREL>>:>)
##############################
###  Enable c++17 support  ###
##############################
target_compile_features(dmrg-flags INTERFACE cxx_std_17)


#######################################
#######################################
# Settings for sanitizers           ###
#######################################
if(DMRG_ENABLE_ASAN)
    target_compile_options(dmrg-flags INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_libraries(dmrg-flags INTERFACE -fsanitize=address)
    if(NOT BUILD_SHARED_LIBS)
        target_link_libraries(dmrg-flags INTERFACE -static-libasan)
    endif()
endif()
if(DMRG_ENABLE_USAN)
    target_compile_options(dmrg-flags INTERFACE -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_libraries(dmrg-flags INTERFACE -fsanitize=undefined)
endif()


#######################################
### Compiler-dependent linker flags ###
#######################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
    target_link_libraries(dmrg-flags INTERFACE  -stdlib=libstdc++ )
endif()


##################################################################
### Check if we can enable Link Time Optimization              ###
##################################################################
if(DMRG_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(DMRG++ PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
        set(DMRG_ENABLE_LTO OFF CACHE)
    endif()
endif()

#################################################################
### Preempt Threads::Threads                                   ###
### It's looked for in dependencies, so we make it right       ###
### before it's done wrong, i.e. with pthread instead of       ###
### -lpthread.                                                 ###
### Here we specify the linking twice                          ###
### 1) As string to make sure -lpthread gets sandwiched by     ###
###    -Wl,--whole-archive.... -Wl,--no-whole-archive          ###
### 2) As usual to make sure that if somebody links            ###
###    Threads::Threads, then any repeated pthread appended    ###
###    to the end (the wrong order causes linking errors)      ###
##################################################################
if(NOT TARGET Threads::Threads)
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    find_package(Threads)
    target_link_libraries(Threads::Threads INTERFACE rt dl)
endif ()


###########################################
### Fetch or find all the dependencies  ###
###########################################
include(cmake/SetupDependencies.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(dmrg-arp PUBLIC dmrg-flags)
target_link_libraries(dmrg-eig PUBLIC dmrg-flags)
target_link_libraries(dmrg-opt PUBLIC dmrg-flags)
target_link_libraries(dmrg-dbg PUBLIC dmrg-flags)
target_link_libraries(dmrg-main PUBLIC dmrg-flags)
target_link_libraries(DMRG++ PUBLIC dmrg-main dmrg-eig dmrg-arp dmrg-opt dmrg-dbg dmrg-deps dmrg-flags)  # Link everything to main project
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "8.0.0")
#    target_compile_options(dmrg-opt PRIVATE -ftime-trace)
#    target_compile_options(dmrg-eig PRIVATE -ftime-trace)
#    target_compile_options(dmrg-arp PRIVATE -ftime-trace)
#    target_compile_options(dmrg-main PRIVATE -ftime-trace)
#    target_compile_options(DMRG++ PRIVATE -ftime-trace)
#endif()
##################################################################
###  Link system libs statically                               ###
##################################################################
if (NOT BUILD_SHARED_LIBS)
    target_link_options(DMRG++ PUBLIC -static-libgcc -static-libstdc++)
endif ()


##################################################################
### Speed up compilation with precompiled headers              ###
##################################################################
if (DMRG_ENABLE_PCH)
    if (COMMAND target_precompile_headers)
        target_precompile_headers(dmrg-flags INTERFACE <string> <vector> <array> <optional> <complex> <memory>)
        target_precompile_headers(dmrg-eig PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor>)
        target_precompile_headers(dmrg-opt PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor>)
        target_precompile_headers(dmrg-main PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor> <chrono> <string_view>)

        target_precompile_headers(dmrg-eig PRIVATE <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
        target_precompile_headers(dmrg-opt PRIVATE <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
        target_precompile_headers(dmrg-main PRIVATE <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
    endif()
endif()

#######################################
# Settings for ccache               ###
#######################################
if(DMRG_ENABLE_CCACHE)
    set(COMPILER_OK ON)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(dmrg-flags INTERFACE -Xclang -fno-pch-timestamp -fpch-preprocess)
    endif()

    mark_as_advanced(COMPILER_OK)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0.0")
        set(COMPILER_OK OFF)
    endif()
    if(COMPILER_OK)
        find_program(CCACHE_PROGRAM ccache)
        if(CCACHE_PROGRAM)
            set_target_properties(DMRG++ PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            set_target_properties(dmrg-main PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            set_target_properties(dmrg-opt PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            set_target_properties(dmrg-dbg PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            set_target_properties(dmrg-arp PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            set_target_properties(dmrg-eig PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
            message(STATUS "Using ccache ${CCACHE_PROGRAM}")
            if(DMRG_ENABLE_PCH AND COMMAND target_precompile_headers)
                message(STATUS "Detected ccache + pch: Remember to set --> sloppiness = include_file_mtime,pch_defines,time_macros <-- in your ccache.conf")
            endif()
        else()
            message(STATUS "ccache program could not be found")
        endif()
    endif()
endif()




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
if (DMRG_PRINT_INFO)
    include(cmake/PrintTargetInfo.cmake)
    include(cmake/getExpandedTarget.cmake)
    message(STATUS "| DMRG TARGET SUMMARY")
    message(STATUS "|--------------------")
    print_compiler_info("| ")
    expand_target_all_targets(DMRG++ DMRG_TARGET_EXPANDED)
    foreach (tgt ${DMRG_TARGET_EXPANDED})
        print_target_info(${tgt} "| ")
    endforeach ()
endif ()


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






