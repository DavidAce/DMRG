cmake_minimum_required(VERSION 3.9)
project(DMRG++)


set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake_modules/gitversion.cmake)


##################################################################
### Define variables for directory paths accessible at runtime ###
##################################################################
include(cmake_modules/directory.cmake)


##################################################################
### Add a debugging type message, only active in Debug-mode   ###
##################################################################
include(cmake_modules/debug_msg.cmake)




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
set(COMMON_OPTIONS -g -Wall -march=native)
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)

##########################
### OS-dependent flags ###
##########################
if(${CMAKE_HOST_APPLE})
    list(APPEND COMMON_OPTIONS -mmacosx-version-min=10.13)
else()
    list(APPEND COMMON_OPTIONS -mavx -mfma)
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
endif()





###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC source)



###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_options (${PROJECT_NAME} PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_link_libraries  (${PROJECT_NAME} PRIVATE -flto)                                               ### Link time optimization flag
target_link_libraries  (${PROJECT_NAME} PRIVATE -lstdc++fs)                                          ### For experimental::filesystem


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

set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)
#include(cmake_modules/Find_dont_install_INTELMKL.cmake)     # MKL - Intel's Math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
include(cmake_modules/Fetch_OpenBLAS.cmake)                 # If MKL is not found openblas will be used instead. Includes lapack.
include(cmake_modules/Fetch_Eigen3.cmake)                   # Eigen library
include(cmake_modules/Fetch_HDF5.cmake)                     # HDF5 binary file format for output data.
include(cmake_modules/Fetch_Arpack-ng.cmake)                # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake_modules/Fetch_Arpackpp.cmake)                 # C++ frontend for arpack-ng
include(cmake_modules/Fetch_GSL.cmake)                      # GSL for computing some numerical integrals.




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






message("=============================================================================")
message("=== SUMMARY OF CMAKE CONFIGURATION:  THE FOLLOWING LIBRARIES WILL BE USED ===")
message("=============================================================================")
message("    INTEL MKL:   ${MKL_LIBRARIES}")
message("    BLAS     :   ${BLAS_LIBRARIES}")
message("    LAPACK   :   ${LAPACK_LIBRARIES}")
message("    EIGEN3   :   ${EIGEN3_INCLUDE_DIR}")
message("    HDF5     :   ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES}")
message("    ARPACK   :   ${ARPACK_LIBRARIES}")
message("    ARPACKPP :   ${ARPACKPP_INCLUDE_DIR}")
message("    GSL      :   ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES}")
message("=============================================================================")
message("")




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
# TEST
############################

enable_testing()

include(tests/arpack++_simple_test.cmake)
include(tests/arpack++_mps_test.cmake)

add_custom_target(all_tests ALL
        DEPENDS arpack++_simple_test_target
        DEPENDS arpack++_mps_test_target
        DEPENDS blas lapack EIGEN3 arpack arpack++ hdf5::hdf5 GSL
        DEPENDS ${PROJECT_NAME}
        )


add_custom_command(
        TARGET all_tests
        COMMENT "Running Tests"
        PRE_BUILD
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failures
        DEPENDS ${PROJECT_NAME}
)


#add_dependencies(${PROJECT_NAME} all_tests)

# To print all variables, use the code below:
#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
