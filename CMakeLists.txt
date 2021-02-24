cmake_minimum_required(VERSION 3.12)
project(DMRG++
        VERSION 2.1.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")


# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake-modules)

# Set options this build
option(DMRG_ENABLE_OPENMP            "Use OpenMP in Eigen, BLAS and LAPACK"                         OFF)
option(DMRG_ENABLE_THREADS           "Use C++11 stl threads in Eigen::Tensor"                       OFF)
option(DMRG_ENABLE_MKL               "Enable Intel Math Kernel Library"                             OFF)
option(DMRG_ENABLE_TESTS             "Enable unit testing with ctest"                               OFF)
option(DMRG_ENABLE_ASAN              "Enable runtime address sanitizer -fsanitize=address"          OFF)
option(DMRG_ENABLE_USAN              "Enable undefined behavior sanitizer -fsanitize=undefined"     OFF)
option(DMRG_ENABLE_LTO               "Enable link time optimization"                                OFF)
option(DMRG_ENABLE_PCH               "Enable precompiled headers to speed up compilation"           ON)
option(DMRG_ENABLE_CCACHE            "Enable ccache to speed up compilation"                        ON)
option(DMRG_BUILD_EXAMPLES           "Build examples"                                               OFF)
option(DMRG_PRINT_INFO               "Print information during cmake configure"                     OFF)
option(DMRG_PRINT_CHECKS             "Print more information during cmake configure"                OFF)
option(DMRG_DEPS_IN_SUBDIR           "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>"     ON)
option(DMRG_PREFER_CONDA_LIBS        "Prefer libraries installed through anaconda"                  OFF)

# Make an "enum" for valid download methods
set(DMRG_PACKAGE_MANAGERS_VALID find cmake find-or-cmake conan)
set(DMRG_PACKAGE_MANAGER find CACHE STRING "Download method for external dependencies")
set_property(CACHE DMRG_PACKAGE_MANAGER PROPERTY STRINGS ${DMRG_PACKAGE_MANAGERS_VALID})
if (NOT DMRG_PACKAGE_MANAGER IN_LIST DMRG_PACKAGE_MANAGERS_VALID)
    message(FATAL_ERROR "DMRG_PACKAGE_MANAGER must be one of ${DMRG_PACKAGE_MANAGERS_VALID}")
endif ()




# Set default policies if CMake is new enough
if (CMAKE_VERSION VERSION_LESS 3.12)
    message(STATUS "Not setting policies")
else()
    cmake_policy(SET CMP0074 NEW) # To let find_package calls use <PackageName>_ROOT variables
    cmake_policy(SET CMP0075 NEW) # Include file check macros honor CMAKE_REQUIRED_LIBRARIES
    # Use cmake_policy CMP0067 to let try_compile honor current CMAKE_CXX_STANDARD setting.
    # See here https://stackoverflow.com/questions/47213356/cmake-using-corrext-c-standard-when-checking-for-header-files
    cmake_policy(SET CMP0067 NEW)
endif()



################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake-modules/gitversion.cmake)


##################################################################
### Setup dependency paths                                     ###
##################################################################
include(GNUInstallDirs)
string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Setting default install prefix -- CMAKE_INSTALL_PREFIX --> ${CMAKE_BINARY_DIR}/dmrg-deps-install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/dmrg-deps-install" CACHE PATH "..." FORCE)
endif()



############################################################
### Print extra info during CMake configure              ###
############################################################
include(cmake-modules/PrintHostInfo.cmake)
include(cmake-modules/PrintOptionInfo.cmake)

###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
include(cmake-modules/CompilerFlags.cmake)



################################
###  Add all source files    ###
################################
include(source_list.txt)
add_executable(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_WHAT_YOU_USE TRUE)
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
        source/math/eig/solution.cpp
        source/math/eig/arpack_solver/matrix_product_dense.cpp
        source/math/eig/arpack_solver/matrix_product_sparse.cpp
        source/math/eig/arpack_solver/matrix_product_hamiltonian.cpp
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
        source/tools/finite/opt-internal/arpack_energy.cpp
        source/tools/finite/opt-internal/arpack_variance.cpp
        source/tools/finite/opt-internal/report.cpp
        source/tools/finite/opt.cpp
        )

add_library(dmrg-dbg OBJECT
        source/general/stack_trace.cpp
        )



target_include_directories(dmrg-main PRIVATE source)
target_include_directories(dmrg-eig  PRIVATE source)
target_include_directories(dmrg-arp  PRIVATE source)
target_include_directories(dmrg-opt  PRIVATE source)
target_include_directories(dmrg-dbg  PRIVATE source)

add_library(dmrg-flags INTERFACE)
target_include_directories(dmrg-flags SYSTEM INTERFACE source)

##################################################################
###  Force cmake to find .a library suffixes on static builds  ###
##################################################################
if(BUILD_SHARED_LIBS)
    set(LINK_TYPE SHARED)
else()
    set(LINK_TYPE STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
    if(NOT DMRG_ENABLE_ASAN)
        target_link_libraries  (${PROJECT_NAME} PUBLIC -static)   ### Static linkage
    endif()
endif()






#######################################################
###  Check  #include<optional> or                   ###
###  or     #include<experimental/optional>         ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckCXXOptional.cmake)
CheckCXXOptional()



##################################################
###  Add optional RELEASE/DEBUG compile flags  ###
##################################################
# -fsanitize=leak -fsanitize=undefined -fsanitize=alignment  -fsanitize=leak -fsanitize=undefined -fsanitize=alignment
target_compile_options(dmrg-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:-g>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>: -Wpedantic -Wextra -Wconversion>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:>)
target_compile_options(dmrg-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:MINSIZEREL>>:>)
##############################
###  Enable c++17 support  ###
##############################
target_compile_features(dmrg-flags INTERFACE cxx_std_17)


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
        set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
        set(DMRG_ENABLE_LTO OFF CACHE)
    endif()
endif()





#################################################################
### Prepend conda dirs for hints when using find_package      ###
#################################################################
if(DMRG_PREFER_CONDA_LIBS)
    list(APPEND CONDA_HINTS
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3
            $ENV{HOME}/miniconda
            )
endif()


##################################################################
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
#    set(THREADS_PREFER_PTHREAD_FLAG FALSE)
    find_package(Threads)
#    set_target_properties(Threads::Threads PROPERTIES INTERFACE_LINK_LIBRARIES "")
    #    target_link_libraries(Threads::Threads INTERFACE "-Wl,--whole-archive;-lpthread;-Wl,--no-whole-archive;rt;dl")
#    target_link_libraries(Threads::Threads INTERFACE -Wl,--whole-archive,-lpthread,--no-whole-archive;rt;dl)
    target_link_libraries(Threads::Threads INTERFACE -Wl,--allow-multiple-definition,--whole-archive,-lpthread,--no-whole-archive;rt;dl)
#        target_link_libraries(Threads::Threads INTERFACE pthread rt dl)
    target_link_libraries(Threads::Threads INTERFACE rt dl)
endif()


###########################################
### Fetch or find all the dependencies  ###
###########################################
include(cmake-modules/SetupSearchPaths.cmake)
include(cmake-modules/SetupDependenciesCMake.cmake)
include(cmake-modules/SetupDependenciesConan.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
if(TARGET openmp::openmp)
    target_link_libraries(dmrg-flags INTERFACE openmp::openmp)
endif()
target_link_libraries(dmrg-arp  PUBLIC dmrg-flags)
target_link_libraries(dmrg-eig  PUBLIC dmrg-flags)
target_link_libraries(dmrg-opt  PUBLIC dmrg-flags)
target_link_libraries(dmrg-dbg  PUBLIC dmrg-flags)
target_link_libraries(dmrg-main PUBLIC dmrg-flags)
target_link_libraries(${PROJECT_NAME} PUBLIC  dmrg-main dmrg-eig dmrg-arp dmrg-opt dmrg-dbg dmrg-flags)  # Link everything to main project
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "8.0.0")
#    target_compile_options(dmrg-opt PRIVATE -ftime-trace)
#    target_compile_options(dmrg-eig PRIVATE -ftime-trace)
#    target_compile_options(dmrg-arp PRIVATE -ftime-trace)
#    target_compile_options(dmrg-main PRIVATE -ftime-trace)
#    target_compile_options(${PROJECT_NAME} PRIVATE -ftime-trace)
#endif()



##################################################################
### Speed up compilation with precompiled headers              ###
##################################################################
if(DMRG_ENABLE_PCH)
    if(COMMAND target_precompile_headers)
        target_precompile_headers(dmrg-flags INTERFACE <string> <vector> <array> <optional> <complex> <memory>)
        target_precompile_headers(dmrg-eig PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor>)
        target_precompile_headers(dmrg-opt PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor>)
        target_precompile_headers(dmrg-main PRIVATE <Eigen/Core> <unsupported/Eigen/CXX11/Tensor> <chrono> <string_view>)

        target_precompile_headers(dmrg-eig PRIVATE  <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
        target_precompile_headers(dmrg-opt PRIVATE  <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
        target_precompile_headers(dmrg-main PRIVATE <config/enums.h> <tools/common/log.h> <tools/common/fmt.h>)
    endif()
endif()

#######################################
# Settings for ccache               ###
#######################################
if(DMRG_ENABLE_CCACHE)
    set(COMPILER_OK ON)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        #        target_compile_options(dmrg-main PRIVATE -Xclang -fno-pch-timestamp)
        #        target_compile_options(dmrg-eig PRIVATE -Xclang -fno-pch-timestamp)
        #        target_compile_options(dmrg-opt PRIVATE -Xclang -fno-pch-timestamp)
    endif()

    mark_as_advanced(COMPILER_OK)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0.0")
        set(COMPILER_OK OFF)
    endif()
    if(COMPILER_OK)
        find_program(CCACHE_PROGRAM ccache)
        if(CCACHE_PROGRAM)
            set_target_properties(${PROJECT_NAME} PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
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
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules/DMRGUninstall.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake
            IMMEDIATE @ONLY)

    add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/DMRGUninstall.cmake)
endif()



# Print summary of CMake configuration
if (DMRG_PRINT_INFO)
    include(cmake-modules/PrintTargetInfo.cmake)
    include(cmake-modules/getExpandedTarget.cmake)
    message(STATUS "| DMRG TARGET SUMMARY")
    message(STATUS "|--------------------")
    print_compiler_info("| ")
    expand_target_all_targets(${PROJECT_NAME} DMRG_TARGET_EXPANDED)
    foreach(tgt ${DMRG_TARGET_EXPANDED})
        print_target_info(${tgt} "| ")
    endforeach()
endif()



###################
### INFORMATION ###
###################
message("")
message("=============================================================================")
message("")
message("      To clean cmake files, run './build.sh -c'")
message("      To see more options,  run './build.sh -h'")
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





