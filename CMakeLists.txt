cmake_minimum_required(VERSION 3.12)
project(DMRG++
        VERSION 2.0.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")


# Make sure we use DMRG's own find modules
list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake-modules)

# Set options this build
option(DMRG_PREFER_CONDA_LIBS        "Prefer libraries installed through anaconda"         OFF)
option(DMRG_ENABLE_OPENMP            "Use OpenMP"                                          OFF)
option(DMRG_ENABLE_MKL               "Enable Intel Math Kernel Library"                    OFF)
option(DMRG_ENABLE_LTO               "Enable link time optimization"                       OFF)
option(DMRG_ENABLE_TESTS             "Enable unit testing with ctest"                      OFF)
option(DMRG_BUILD_EXAMPLES           "Build examples"                                      OFF)
option(DMRG_PRINT_INFO               "Print information during cmake configure"            OFF)



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
    cmake_policy(SET CMP0074 NEW)
    cmake_policy(SET CMP0075 NEW)
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



############################################################
### Set  the same microarchitecture for c++ and OpenBLAS ###
############################################################

if(NOT DMRG_MARCH)
    set(DMRG_MARCH "native")
endif()
if(DMRG_MARCH)
    if (${DMRG_MARCH} STREQUAL "zen")
        string(TOUPPER ${DMRG_MARCH} OPENBLAS_MARCH)
        set(CXX_MARCH znver1)
    elseif (${DMRG_MARCH} STREQUAL "native")
        set(OPENBLAS_MARCH HASWELL)
        set(CXX_MARCH native)
    else()
        string(TOUPPER ${DMRG_MARCH} OPENBLAS_MARCH)
        string(TOLOWER ${DMRG_MARCH} CXX_MARCH)
    endif()
endif()




################################
###  Add all source files    ###
################################
include(source_list.txt)
target_include_directories(${PROJECT_NAME} PUBLIC source)


add_library(project-settings INTERFACE)
##################################################################
###  Force cmake to find .a library suffixes on static builds  ###
##################################################################
if(BUILD_SHARED_LIBS)
    set(LINK_TYPE SHARED)
else()
    set(LINK_TYPE STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
    target_link_libraries  (project-settings INTERFACE -static)   ### Static linkage
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



###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
# Segfault from -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC when using ceres, but -D_GLIBCXX_DEBUG fixes gdb incomplete types
message(STATUS "Using microarchitechture: ${CXX_MARCH}")
target_compile_options(project-settings INTERFACE $<$<COMPILE_LANGUAGE:CXX>: -Wall -march=${CXX_MARCH} -mtune=${CXX_MARCH}>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:          -O3 -g -DNDEBUG>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:            -O0 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:   -O1 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:PROFILE>>:          -O3 -g  -DNDEBUG  -lprofiler  -ftime-report>)
target_compile_features(project-settings INTERFACE cxx_std_17)

#######################################
### Compiler-dependent linker flags ###
#######################################
target_link_libraries(project-settings INTERFACE  -lstdc++fs)
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
include(cmake-modules/SetupDependenciesFind.cmake)
include(cmake-modules/SetupDependenciesFetch.cmake)
include(cmake-modules/SetupDependenciesConan.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(project-settings INTERFACE  ${CONAN_TARGETS} ${NATIVE_TARGETS}) # Order is important here


if(NOT TARGET openmp::openmp)
    target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
endif()
target_link_libraries(project-settings INTERFACE -Wl,--whole-archive pthread -Wl,--no-whole-archive -lrt -ldl )
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)  # Link everything to main project






# Print summary of CMake configuration
if (DMRG_PRINT_INFO)
    include(cmake-modules/PrintTargetInfo.cmake)
    include(cmake-modules/getExpandedTarget.cmake)
    expand_target_all_targets(project-settings DMRG_TARGETS_EXPANDED)
    message(STATUS "| DMRG TARGET SUMMARY")
    message(STATUS "|--------------------")
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
#
#if(DMRG_ENABLE_TESTS)
#    enable_testing()
#    add_subdirectory(tests/arpack++_simple)
#    add_subdirectory(tests/arpack++_mps   )
#    add_subdirectory(tests/file_io        )
#    add_subdirectory(tests/eigsolver      )
#    add_subdirectory(tests/svdwrapper     )
#    add_subdirectory(tests/eigen_openmp   )
#    add_subdirectory(tests/ceres_pthread  )
#
#    add_custom_target(all-tests
#            DEPENDS test-arpack++_simple
#            DEPENDS test-arpack++_mps
#            DEPENDS test-file_io
#            DEPENDS test-eigsolver
#            DEPENDS test-eigen_openmp
#            DEPENDS test-svdwrapper
#            DEPENDS test-ceres_pthread
#            )
#
#    add_custom_command(TARGET all-tests
#            POST_BUILD
#            COMMENT "Running Tests"
#            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --output-on-failure)
#endif()





