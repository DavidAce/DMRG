cmake_minimum_required(VERSION 3.10)
project(DMRG++)


set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake-modules)


################################################################
### Set options for compilation                              ###
################################################################

#option(STATIC_BUILD  "Static linking"                   ON)

################################################################
### Set options for compilation                              ###
################################################################
option(BUILD_SHARED_LIBS  "Sets shared linking on/off"                          OFF)
option(BUILD_EXAMPLES     "Build examples"                                      OFF)
option(ENABLE_TESTS       "Enable testing"                                      ON)
option(DOWNLOAD_HDF5      "Download and install HDF5 if not found"              OFF)
option(DOWNLOAD_SPDLOG    "Download and install spdlog if not found"            OFF)
option(DOWNLOAD_EIGEN3    "Download and install Eigen3 if not found"            OFF)
option(DOWNLOAD_ALL       "Download and install all dependencies if not found"  ON)
option(ENABLE_CCACHE      "Use ccache to speed up compilation"                  OFF)
option(ENABLE_COTIRE      "Use cotire to speed up compilation"                  OFF)
option(USE_OpenMP         "Use OpenMP"                                          OFF)
option(USE_MKL            "Use Intel MKL libraries"                             ON)



#################################################################
### Set default policies if CMake is new enough               ###
#################################################################
if (CMAKE_VERSION VERSION_LESS 3.12)
    message(STATUS "Not setting policies")
else()
    cmake_policy(SET CMP0074 NEW)
    cmake_policy(SET CMP0075 NEW)
endif()


if(USE_MKL AND NOT DEFINED ENV{MKLROOT})
    set(ENV{MKLROOT}  "/opt/intel/mkl")
    message(STATUS "MKLROOT is not defined. Setting default: $ENV{MKLROOT}")
endif()


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake-modules/gitversion.cmake)


##################################################################
### Define variables for directory paths accessible at runtime ###
##################################################################
#include(cmake-modules/directory.cmake)


##################################################################
### Add a debugging type message, only active in Debug-mode   ###
##################################################################
#include(cmake-modules/debug_msg.cmake)




##################################################################
### Get operating system properties into a variable            ###
##################################################################
set(GET_OS_INFO_CMD lsb_release -a)
if(${CMAKE_HOST_APPLE})
    set(GET_OS_INFO_CMD "sw_vers")
endif()
execute_process(COMMAND ${GET_OS_INFO_CMD}
        OUTPUT_VARIABLE OS_PROPERTIES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message("========= DETECTED OS =========")
message("${OS_PROPERTIES}")
message("===============================")


###################
### Setup flags ###
###################
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if(NOT MARCH)
    set(MARCH "native")
endif()

if(MARCH)
    if (${MARCH} STREQUAL "zen")
        string(TOUPPER ZEN OPENBLAS_MARCH)
        string(TOLOWER znver1 CXX_MARCH)
    elseif (${MARCH} STREQUAL "native")
        set(OPENBLAS_MARCH "haswell")
        string(TOLOWER ${MARCH} CXX_MARCH)
    else()
        string(TOUPPER ${MARCH} OPENBLAS_MARCH)
        string(TOLOWER ${MARCH} CXX_MARCH)
    endif()
endif()

message("Using microarchitechture: ${CXX_MARCH}")
if(ENABLE_CCACHE)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
        set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
    endif(CCACHE_FOUND)
endif()
set(COMMON_OPTIONS -g -Wall -march=${CXX_MARCH} -mtune=${CXX_MARCH})
set(DEBUG_OPTIONS       -Wpedantic -Wextra -O0 )
set(RELEASE_OPTIONS     -DNDEBUG -O3)
set(PROFILE_OPTIONS     -DNDEBUG -O3 -lprofiler -g  -ftime-report)

#################################
### Microarch dependent flags ###
#################################
# Check supported options here http://man7.org/linux/man-pages/man1/gcc.1.html
#if (${MARCH} STREQUAL "sandybridge")
#    list(APPEND RELEASE_OPTIONS  -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -msse4.1  -mpopcnt  -maes -mpclmul -mavx -fstack-protector-strong -Wformat -Wformat-security)
#elseif(${MARCH} STREQUAL "haswell" OR ${MARCH} STREQUAL "skylake" OR ${MARCH} STREQUAL "zen" )
#    list(APPEND RELEASE_OPTIONS -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -msse4.1 -mpopcnt -maes -mpclmul -mavx -mavx2 -mmovbe -mrdrnd -mf16c -mfsgsbase -mfma  -mbmi -mbmi2 -fstack-protector-strong -Wformat -Wformat-security)
#else()
#    # Set default to whatever nehalem can manage
#    list(APPEND RELEASE_OPTIONS  -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -msse4.1 -mpopcnt -fstack-protector-strong -Wformat -Wformat-security)
#endif()


##########################
### OS-dependent flags ###
##########################
if(${CMAKE_HOST_APPLE})
    list(APPEND COMMON_OPTIONS -mmacosx-version-min=10.13)
else()
#    list(APPEND COMMON_OPTIONS -mavx ) #-mfma makes things slower?
endif()


################################
### Compiler-dependent flags ###
################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    list(APPEND COMMON_OPTIONS -Wno-implicit-fallthrough -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-enum-compare)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    list(APPEND COMMON_OPTIONS -stdlib=libstdc++ -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-invalid-partial-specialization -Wno-missing-braces -Wno-overloaded-virtual -Wno-uninitialized)
endif()


################################
### Print current CXX flags  ###
################################
message(STATUS "COMMON CXX flags:  ${COMMON_OPTIONS} ${CMAKE_CXX_FLAGS}")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "RELEASE CXX flags: ${RELEASE_OPTIONS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "DEBUG CXX flags: ${DEBUG_OPTIONS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
    message(STATUS "PROFILE CXX flags: ${PROFILE_OPTIONS}")
endif()





###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC source)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(BUILD_SHARED_LIBS)
    set(CUSTOM_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
else()
    set(CUSTOM_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
    target_link_libraries  (${PROJECT_NAME} PRIVATE -static)                                             ### Static linkage
endif()



#endif()



##############################################
###  Optional OpenMP support               ###
###  Note that Clang has some              ###
###  trouble with static openmp            ###
###  Also, installing libiomp5 might help  ###
###  for shared linking.                   ###
##############################################
if(USE_OpenMP)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND STATIC_BUILD)
        message(WARNING "Disabling OpenMP
        Clang cannot link \"libomp\" statically.
        Please use GNU C++ instead to link \"libgomp.a\" statically, or turn off OpenMP
        by passing \"-DUSE_OpenMP:BOOL=OFF\" to CMake to avoid this warning.")
        message(STATUS "OpenMP OFF")
        set(USE_OpenMP OFF)
    else()
        find_package(OpenMP)
        if (OpenMP_FOUND)
            message(STATUS "OpenMP ON" )
            set(OpenMP_LIBRARIES  ${OpenMP_gomp_LIBRARY} ${OpenMP_omp_LIBRARY} ${OpenMP_iomp_LIBRARY})
            target_compile_options    (${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
            target_link_libraries     (${PROJECT_NAME} PRIVATE ${OpenMP_LIBRARIES})
            target_compile_definitions(${PROJECT_NAME} PRIVATE OpenMP_AVAILABLE)
        else()
            message(WARNING "OpenMP OFF -- Could not be found." )
        endif()
    endif()
#
#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        if("${_variableName}" MATCHES "openmp" OR "${_variableName}" MATCHES "OpenMP" OR "${_variableName}" MATCHES "omp")
#            message(STATUS "${_variableName}=${${_variableName}}")
#        endif()
#    endforeach()
else()
    list(APPEND COMMON_OPTIONS -Wno-unknown-pragmas)
endif()



###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_options (${PROJECT_NAME} PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:PROFILE>:${PROFILE_OPTIONS}>")           ### Profiling build options
target_link_libraries  (${PROJECT_NAME} PRIVATE -flto)                                               ### Link time optimization flag
target_link_libraries  (${PROJECT_NAME} PRIVATE -lstdc++fs)                                          ### For experimental::filesystem
set(PTHREAD_LIBRARY  -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
#target_link_libraries  (${PROJECT_NAME} PRIVATE ${PTHREAD_LIBRARY} )
#target_link_libraries  (${PROJECT_NAME} PRIVATE -Wl,--no-allow-shlib-undefined )
#target_link_libraries  (${PROJECT_NAME} PRIVATE -Wl,--no-undefined )


###############################
###  Require c++17 support  ###
###############################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME}
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)



###########################################
###  Find or fetch all dependencies   ###
###########################################
set(BUILD_DIRECTORY   ${PROJECT_SOURCE_DIR}/cmake-build-libs)
set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)
#include(cmake-modules/Find_dont_install_INTELMKL.cmake)    # MKL - Intel's Math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
include(cmake-modules/FindGFortran.cmake)
include(cmake-modules/Fetch_OpenBLAS.cmake)                 # If MKL is not found openblas will be used instead. Includes lapack.
include(cmake-modules/Fetch_Eigen3.cmake)                   # Eigen library
include(cmake-modules/Fetch_spdlog.cmake)                   # spdlog library
include(cmake-modules/Fetch_LBFGSpp.cmake)                  # Limited memory BFGS gradient descent library
include(cmake-modules/Fetch_Arpack-ng.cmake)                # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake-modules/Fetch_Arpackpp.cmake)                 # C++ frontend for arpack-ng
include(cmake-modules/Fetch_GSL.cmake)                      # GSL for computing some numerical integrals.
include(cmake-modules/Fetch_HDF5.cmake)                     # HDF5 binary file format for output data.
include(cmake-modules/Fetch_h5pp.cmake)                     # GSL for computing some numerical integrals.

target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog)
target_link_libraries(${PROJECT_NAME} PRIVATE LBFGSpp)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5)
target_link_libraries(${PROJECT_NAME} PRIVATE arpack++)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL)
target_link_libraries(${PROJECT_NAME} PRIVATE h5pp::h5pp)

######################################
###  Require experimental headers  ###
######################################
# Use cmake_policy CMP0067 to let try_compile honor current CMAKE_CXX_STANDARD setting.
# See here
# https://stackoverflow.com/questions/47213356/cmake-using-corrext-c-standard-when-checking-for-header-files
cmake_policy(SET CMP0067 NEW)
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_FLAGS -std=c++17  -lstdc++fs -l${PATH_TO_LIBSTDCXXFS})
check_include_file_cxx(experimental/filesystem    has_experimental_filesystem)
check_include_file_cxx(experimental/type_traits   has_experimental_type_traits)
if(NOT has_experimental_filesystem OR NOT has_experimental_type_traits )
    message(FATAL_ERROR "\n\
        Missing one or more experimental headers.\n\
        Consider using a newer compiler (GCC 7 or above, Clang 6 or above),\n\
        or checking the compiler flags.")
endif()



get_target_property(BLAS_INFO       blas             INTERFACE_LINK_LIBRARIES)
get_target_property(LAPACK_INFO     lapack           INTERFACE_LINK_LIBRARIES)
get_target_property(LAPACKE_INFO    lapacke          INTERFACE_LINK_LIBRARIES)
get_target_property(ARPACK_INFO     arpack           INTERFACE_LINK_LIBRARIES)
get_target_property(ARPACKPP_INFO   arpack++         INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(EIGEN_INFO      Eigen3::Eigen    INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LBFGS_INFO      LBFGSpp          INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(SPDLOG_INFO     spdlog::spdlog   INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(HDF5_INFO       hdf5             INTERFACE_LINK_LIBRARIES)
get_target_property(GSL_INFO        GSL              INTERFACE_LINK_LIBRARIES)
get_target_property(FORTRAN_INFO    gfortran         INTERFACE_LINK_LIBRARIES)
get_target_property(H5PP_INFO       h5pp::h5pp       INTERFACE_INCLUDE_DIRECTORIES)



message("=============================================================================")
message("=== SUMMARY OF CMAKE CONFIGURATION:  THE FOLLOWING LIBRARIES WILL BE USED ===")
message("=============================================================================")
message("    BLAS     :   ${BLAS_INFO}")
message("    LAPACK   :   ${LAPACK_INFO}")
message("    LAPACKE  :   ${LAPACKE_INCLUDE_DIRS} ${LAPACKE_INFO}")
message("    EIGEN3   :   ${EIGEN_INFO}")
message("    LBFGS++  :   ${LBFGS_INFO}")
message("    SPDLOG   :   ${SPDLOG_INFO}")
message("    HDF5     :   ${HDF5_INFO}")
message("    ARPACK   :   ${ARPACK_INFO}")
message("    ARPACKPP :   ${ARPACKPP_INFO}")
message("    GSL      :   ${GSL_INFO}")
message("    h5pp     :   ${H5PP_INFO}")
message("    GFORTRAN :   ${FORTRAN_INFO}")
#message("    QUADMATH :   ${QUADMATH_LIB}")
message("=============================================================================")
message("")



###########################################
### Use precompiled headers with cotire ###
###########################################
if (ENABLE_COTIRE)
    get_target_property(LAPACKE_DIRS   lapacke         INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(ARPACK++_DIRS  arpack++        INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COTIRE_PREFIX_HEADER_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/source"
            COTIRE_PREFIX_HEADER_IGNORE_PATH
            "${LAPACKE_DIRS};${ARPACK++_DIRS}"
            )
    include(cmake-modules/cotire.cmake)
    cotire(${PROJECT_NAME})

    #    set_property(DIRECTORY
#            PROPERTY
#            COTIRE_PREFIX_HEADER_IGNORE_PATH
#            ${PROJECT_SOURCE_DIR}/libs/arpack++
#            /usr/include/lapacke.h
#            )


endif()



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


############################
# UNIT TEST
############################

enable_testing()

include(tests/arpack++_simple_test.cmake)
include(tests/arpack++_mps_test.cmake)
include(tests/hdf5_test.cmake)
include(tests/lapacke_test.cmake)

add_custom_target(all_tests ALL
        DEPENDS arpack++_simple_test_target
        DEPENDS arpack++_mps_test_target
        DEPENDS hdf5_test_target
        DEPENDS lapacke_test_target
        DEPENDS arpack++ h5pp::h5pp hdf5 Eigen3::Eigen spdlog::spdlog
        DEPENDS ${PROJECT_NAME}
        )


add_custom_command(
        TARGET all_tests
        COMMENT "Running Tests"
        PRE_BUILD
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failures
        DEPENDS ${PROJECT_NAME} all_tests
)


#add_dependencies(${PROJECT_NAME} all_tests)

# To print all variables, use the code below:
#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
