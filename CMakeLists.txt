cmake_minimum_required(VERSION 3.12)
project(DMRG++
        VERSION 2.1.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")


# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake-modules)

# Set options this build
option(DMRG_ENABLE_OPENMP            "Use OpenMP"                                                   ON)
option(DMRG_ENABLE_MKL               "Enable Intel Math Kernel Library"                             OFF)
option(DMRG_ENABLE_TESTS             "Enable unit testing with ctest"                               OFF)
option(DMRG_ENABLE_ASAN              "Enable runtime address sanitizer -fsanitize=address"          OFF)
option(DMRG_ENABLE_USAN              "Enable undefined behavior sanitizer -fsanitize=undefined"     OFF)
option(DMRG_ENABLE_LTO               "Enable link time optimization"                                OFF)
option(DMRG_ENABLE_PCH               "Enable precompiled headers to speed up compilation"           ON)
option(DMRG_ENABLE_CCACHE            "Enable ccache to speed up compilation"                        ON)
option(DMRG_BUILD_EXAMPLES           "Build examples"                                               OFF)
option(DMRG_PRINT_INFO               "Print information during cmake configure"                     OFF)
option(DMRG_DEPS_IN_SUBDIR           "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>"     ON)
option(DMRG_PREFER_CONDA_LIBS        "Prefer libraries installed through anaconda"                  OFF)

# Make an "enum" for valid download methods
set(DMRG_DOWNLOAD_METHODS_VALID find fetch find-or-fetch conan)
set(DMRG_DOWNLOAD_METHOD find CACHE STRING "Download method for external dependencies")
set_property(CACHE DMRG_DOWNLOAD_METHOD PROPERTY STRINGS ${DMRG_DOWNLOAD_METHODS_VALID})
if (NOT DMRG_DOWNLOAD_METHOD IN_LIST DMRG_DOWNLOAD_METHODS_VALID)
    message(FATAL_ERROR "DMRG_DOWNLOAD_METHOD must be one of ${DMRG_DOWNLOAD_METHODS_VALID}")
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
add_executable(${PROJECT_NAME} ${DMRG_SOURCES} ${DMRG_HEADERS})
add_library(dmrg-eig OBJECT
        source/math/eig/settings.cpp
        source/math/eig/solver.cpp
        source/math/eig/solver_dgeev.cpp
        source/math/eig/solver_dsyevd.cpp
        source/math/eig/solver_zgeev.cpp
        source/math/eig/solver_zheevd.cpp
        source/math/eig/solution.cpp
        source/math/eig/arpack_solver/arpack_solver.cpp
        source/math/eig/arpack_solver/matrix_product_dense.cpp
        source/math/eig/arpack_solver/matrix_product_sparse.cpp
        source/math/eig/arpack_solver/matrix_product_hamiltonian.cpp

        source/math/eig.h
        source/math/eig/settings.h
        source/math/eig/solver.h
        source/math/eig/solution.h
        source/math/eig/view.h
        source/math/eig/arpack_solver/arpack_solver.h
        source/math/eig/arpack_solver/matrix_product_dense.h
        source/math/eig/arpack_solver/matrix_product_sparse.h
        source/math/eig/arpack_solver/matrix_product_hamiltonian.h
        )

target_include_directories(${PROJECT_NAME} PRIVATE source)
target_include_directories(dmrg-eig PRIVATE source)


add_library(project-settings INTERFACE)
##################################################################
###  Force cmake to find .a library suffixes on static builds  ###
##################################################################
if(BUILD_SHARED_LIBS)
    set(LINK_TYPE SHARED)
else()
    set(LINK_TYPE STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
    if(NOT DMRG_ENABLE_ASAN)
        target_link_libraries  (project-settings INTERFACE -static)   ### Static linkage
    endif()
endif()




##############################
###  Enable c++17 support  ###
##############################
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


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
target_compile_options(project-settings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:-g>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>: -Wpedantic -Wextra -Wconversion>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:MINSIZEREL>>:>)
target_compile_features(project-settings INTERFACE cxx_std_17)


#######################################
# Settings for sanitizers           ###
#######################################
if(DMRG_ENABLE_ASAN)
    target_compile_options(project-settings INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_libraries(project-settings INTERFACE -fsanitize=address)
    if(NOT BUILD_SHARED_LIBS)
        target_link_libraries(project-settings INTERFACE -static-libasan)
    endif()
endif()
if(DMRG_ENABLE_USAN)
    target_compile_options(project-settings INTERFACE -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_libraries(project-settings INTERFACE -fsanitize=undefined)
endif()


#######################################
### Compiler-dependent linker flags ###
#######################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
    target_link_libraries(project-settings INTERFACE  -stdlib=libstdc++ )
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
### -lpthread.                                     ###
##################################################################
if(NOT TARGET Threads::Threads)
    set(CMAKE_THREAD_PREFER_PTHREAD FALSE)
    set(THREADS_PREFER_PTHREAD_FLAG FALSE)
    find_package(Threads)
    set_target_properties(Threads::Threads PROPERTIES INTERFACE_LINK_LIBRARIES pthread)
endif()



###########################################
### Fetch or find all the dependencies  ###
###########################################
include(cmake-modules/SetupSearchPaths.cmake)
include(cmake-modules/SetupDependenciesFetch.cmake)
include(cmake-modules/SetupDependenciesConan.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
find_package(Filesystem COMPONENTS Final Experimental REQUIRED)
target_link_libraries(project-settings INTERFACE  std::filesystem)
target_link_libraries(project-settings INTERFACE  ${CONAN_TARGETS} ${FOUND_TARGETS}) # Order is important here

if(NOT TARGET openmp::openmp)
    target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
endif()



# Speed up with precompiled headers
if(DMRG_ENABLE_PCH)
    if(COMMAND target_precompile_headers)
        target_precompile_headers(${PROJECT_NAME} PRIVATE <string> <vector> <array> <list> <set> <complex> <map> <optional> <memory>)
        target_precompile_headers(${PROJECT_NAME} PRIVATE <h5pp/h5pp.h>)
    endif()
endif()

target_link_libraries(project-settings INTERFACE -Wl,--whole-archive pthread -Wl,--no-whole-archive -lrt -ldl )
target_link_libraries(dmrg-eig PRIVATE project-settings )  # Link everything to main project
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-eig project-settings )  # Link everything to main project


#######################################
# Settings for ccache               ###
#######################################
if(DMRG_ENABLE_CCACHE)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        set_target_properties(${PROJECT_NAME} PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        set_target_properties(dmrg-eig PROPERTIES CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        message(STATUS "Using ccache ${CCACHE_PROGRAM}")
        if(DMRG_ENABLE_PCH AND COMMAND target_precompile_headers)
            message(STATUS "Detected ccache + pch: Remember to set --> sloppiness = include_file_mtime,pch_defines,time_macros <-- in your ccache.conf")
        endif()
    else()
        message(STATUS "ccache program could not be found")
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
    expand_target_all_targets(project-settings DMRG_TARGETS_EXPANDED)
    message(STATUS "| DMRG TARGET SUMMARY")
    message(STATUS "|--------------------")
    print_compiler_info("| ")
    foreach(tgt ${DMRG_TARGETS_EXPANDED})
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





