cmake_minimum_required(VERSION 3.12)
project(DMRG++
        VERSION 2.0.0
        DESCRIPTION "DMRG algorithms for 1D quantum spin chains"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake-modules)


################################################################
### Set options for compilation                              ###
################################################################
option(PREFER_CONDA_LIBS        "Prefer libraries installed through anaconda"         OFF)
option(ENABLE_OPENMP            "Use OpenMP"                                          OFF)
option(ENABLE_MKL               "Enable Intel Math Kernel Library"                    OFF)
option(ENABLE_LTO               "Enable link time optimization"                       OFF)
option(BUILD_SHARED_LIBS        "Sets shared linking on/off"                          OFF)
option(BUILD_TESTS              "Builds ctests"                                       OFF)
option(BUILD_EXAMPLES           "Build examples"                                      OFF)
option(DMRG_PRINT_INFO          "Print information during cmake configure"            OFF)


#################################################################
### Make an "enum" for valid download methods                 ###
#################################################################
set(VALID_DOWNLOAD_METHODS none find-only conan native)
set(DOWNLOAD_METHOD find-only CACHE STRING "Download method for external dependencies")
set_property(CACHE DOWNLOAD_METHOD  PROPERTY STRINGS ${VALID_DOWNLOAD_METHODS})
if(NOT DOWNLOAD_METHOD IN_LIST VALID_DOWNLOAD_METHODS)
    message(FATAL_ERROR "DOWNLOAD_METHOD must be one of ${VALID_DOWNLOAD_METHODS}")
endif()

#################################################################
### Set default policies if CMake is new enough               ###
#################################################################
if (CMAKE_VERSION VERSION_LESS 3.12)
    message(STATUS "Not setting policies")
else()
    cmake_policy(SET CMP0074 NEW)
    cmake_policy(SET CMP0075 NEW)
    # Use cmake_policy CMP0067 to let try_compile honor current CMAKE_CXX_STANDARD setting.
    # See here https://stackoverflow.com/questions/47213356/cmake-using-corrext-c-standard-when-checking-for-header-files
    cmake_policy(SET CMP0067 NEW)
endif()

# Set MKLROOT if not defined in environment
if(ENABLE_MKL AND NOT MKLROOT)
    if(DEFINED $ENV{MKLROOT})
        set(MKLROOT $ENV{MKLROOT})
    else()
    set(MKLROOT  /opt/intel/mkl)
    message(STATUS "MKLROOT is not defined. Setting default: ${MKLROOT}")
    endif()
endif()



################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake-modules/gitversion.cmake)




##################################################################
### Print host properties                                      ###
##################################################################
cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
cmake_host_system_information(RESULT _proc_type   QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT _os_name     QUERY OS_NAME)
cmake_host_system_information(RESULT _os_release  QUERY OS_RELEASE)
cmake_host_system_information(RESULT _os_version  QUERY OS_VERSION)
cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)
message("========= HOST INFO =========")
message("${_host_name}")
message("${_os_name} ${_os_platform} ${_os_release}")
message("${_proc_type}")
message("${_os_version}")
message("===============================")



###########################
### Print compiler info ###
###########################
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")


############################################################
### Set  the same microarchitecture for c++ and OpenBLAS ###
############################################################

if(NOT MARCH)
    set(MARCH "native")
endif()
if(MARCH)
    if (${MARCH} STREQUAL "zen")
        string(TOUPPER ${MARCH} OPENBLAS_MARCH)
        set(CXX_MARCH znver1)
    elseif (${MARCH} STREQUAL "native")
        set(OPENBLAS_MARCH HASWELL)
        set(CXX_MARCH native)
    else()
        string(TOUPPER ${MARCH} OPENBLAS_MARCH)
        string(TOLOWER ${MARCH} CXX_MARCH)
    endif()
endif()




################################
###  Add all source files    ###
################################
include(source_list.txt)
target_include_directories(${PROJECT_NAME} PUBLIC source)


add_library(project-settings INTERFACE) # A helper target

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
###  Check  #include<filesystem> or                 ###
###  or     #include<experimental/filesystem>       ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckCXXFilesystem.cmake)
CheckCXXFilesystem()

#######################################################
###  Check  #include<experimental/type_traits>      ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckTypeTraits.cmake)
CheckTypeTraits()

#######################################################
###  Check  #include<optional> or                   ###
###  or     #include<experimental/optional>         ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckCXXOptional.cmake)
CheckCXXOptional()



###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
message(STATUS "Using microarchitechture: ${CXX_MARCH}")
target_compile_options(project-settings INTERFACE $<$<COMPILE_LANGUAGE:CXX>: -Wall -march=${CXX_MARCH} -mtune=${CXX_MARCH}>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:          -O3 -g -DNDEBUG>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:            -O0 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2>) # Segfault from -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC when using ceres
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:   -O1 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2>) # Segfault from -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC when using ceres
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
if(ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
        set(ENABLE_LTO OFF CACHE)
    endif()
endif()


##################################################################
### Setup dependency paths                                     ###
##################################################################
include(GNUInstallDirs)
string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Setting default install prefix -- CMAKE_INSTALL_PREFIX --> ${CMAKE_BINARY_DIR}/dmrg-deps-install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/dmrg-deps-install" CACHE PATH "..." FORCE)
endif()
#set(EXTERNAL_INSTALL_DIR       ${CMAKE_INSTALL_PREFIX})
#set(EXTERNAL_BUILD_DIR         ${CMAKE_BINARY_DIR}/dmrg-deps-build)
#set(EXTERNAL_CONFIG_DIR        ${CMAKE_BINARY_DIR}/dmrg-deps-build)




#################################################################
### Prepend conda dirs for hints when using find_package      ###
#################################################################
if(PREFER_CONDA_LIBS)
    list(APPEND CONDA_HINTS
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3
            $ENV{HOME}/miniconda
            )
endif()

##################################################################
### Print summary of CMake configuration                       ###
##################################################################
message("=============================================================================")
message("===               DMRG++: SUMMARY OF CMAKE CONFIGURATION                  ===")
message("=============================================================================")
message(STATUS "CMAKE_INSTALL_PREFIX         : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "PREFER_CONDA_LIBS            : ${PREFER_CONDA_LIBS}")
message(STATUS "ENABLE_OPENMP                : ${ENABLE_OPENMP}")
message(STATUS "ENABLE_MKL                   : ${ENABLE_MKL}")
message(STATUS "ENABLE_LTO                   : ${ENABLE_LTO}")
message(STATUS "BUILD_SHARED_LIBS            : ${BUILD_SHARED_LIBS}")
message(STATUS "BUILD_TESTS                  : ${BUILD_TESTS}")
message(STATUS "BUILD_EXAMPLES               : ${BUILD_EXAMPLES}")
message("=============================================================================")


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
if(DOWNLOAD_METHOD MATCHES "conan")
    include(cmake-modules/SetupDependenciesConan.cmake)
elseif(DOWNLOAD_METHOD MATCHES "find|native")
    include(cmake-modules/SetupDependenciesNative.cmake)
endif()
include(cmake-modules/PrintTargetInfo.cmake)
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)  # Link everything to main project
print_target_info(${PROJECT_NAME})

###################
### INFORMATION ###
###################
message("")
message("=============================================================================")
message("")
message("      To clean cmake files, run './build.sh -c'")
message("      To clean cmake files downloaded libraries, run './build.sh -l '")
message("      To launch this program,  './run.sh'")
message("")
message("=============================================================================")
message("")



#############################
## UNIT TEST
#############################
if(ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests/arpack++_simple)
    add_subdirectory(tests/arpack++_mps   )
    add_subdirectory(tests/file_io        )
    add_subdirectory(tests/eigsolver      )
    add_subdirectory(tests/svdwrapper     )
    add_subdirectory(tests/eigen_openmp   )
    add_subdirectory(tests/ceres_pthread  )

    add_custom_target(all-tests
            DEPENDS test-arpack++_simple
            DEPENDS test-arpack++_mps
            DEPENDS test-file_io
            DEPENDS test-eigsolver
            DEPENDS test-eigen_openmp
            DEPENDS test-svdwrapper
            DEPENDS test-ceres_pthread
            )

    add_custom_command(TARGET all-tests
            POST_BUILD
            COMMENT "Running Tests"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failure)
endif()





