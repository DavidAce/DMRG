cmake_minimum_required(VERSION 3.8)
project(DMRG++)


set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(cmake/get_libraries.cmake)

################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake/gitversion.cmake)
include(cmake/directory.cmake)

message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")


###################
### Setup flags ###
###################
set(COMMON_OPTIONS -g -Wall -march=native)
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    set(COMMON_OPTIONS ${COMMON_OPTIONS}  -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-maybe-uninitialized)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(COMMON_OPTIONS ${COMMON_OPTIONS}   -stdlib=libstdc++ -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-invalid-partial-specialization -Wno-missing-braces -Wno-overloaded-virtual -Wno-uninitialized)
endif()
message(STATUS "Default CXX flags:  ${COMMON_OPTIONS} ${CMAKE_CXX_FLAGS}")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Additional Release flags: ${RELEASE_OPTIONS}")
else()
    message(STATUS "Additional Debug flags: ${DEBUG_OPTIONS}")
endif()


###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_options (${PROJECT_NAME} PUBLIC ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (${PROJECT_NAME} PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_link_libraries  (${PROJECT_NAME} -flto)                                                       ### Link time optimization flag

######################################################
###  Require c++17 support and test that it works  ###
######################################################
set_target_properties  (${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)                                         ### Demand c++17 compiler support
INCLUDE(CheckIncludeFileCXX)
check_include_file_cxx(experimental/filesystem  has_experimental_filesystem -std=c++17)
check_include_file_cxx(experimental/type_traits has_experimental_type_traits -std=c++17)
if(has_experimental_filesystem AND has_experimental_type_traits)
    message(STATUS "C++17 experimental features are supported.")
else()
    message(FATAL_ERROR "Missing experimental headers for C++17. Check your compiler and libstdc++ versions")
endif()




#######################################
###  Link and include dependencies  ###
#######################################

target_include_directories(${PROJECT_NAME} PRIVATE source)
target_link_libraries(${PROJECT_NAME} -lstdc++fs)                              ### For experimental::filesystem


###########################################
###  Find or install all dependencies   ###
###########################################
#include(cmake/Find_dont_install_INTELMKL.cmake)          ### For Intel MKL library. Note that it conflicts with arpack++
include(cmake/Find_or_install_EIGEN3.cmake)              ### For Eigen linear algebra header-library http://eigen.tuxfamily.org/
include(cmake/Find_or_install_HDF5.cmake)                ### For HDF5 data env_storage format            https://support.hdfgroup.org
include(cmake/Find_or_install_LAPACK_and_BLAS.cmake)     ### For Lapack and Blas libraries           https://github.com/Reference-LAPACK/lapack
include(cmake/Find_or_install_ARPACK.cmake)              ### For Arpack-ng, a mantained arpack fork  https://github.com/opencollab/arpack-ng
include(cmake/Find_or_install_ARPACKPP.cmake)            ### For ARPACKPP header-arplibrary          https://github.com/m-reuter/arpackpp
include(cmake/Find_or_install_GSL.cmake)                 ### For GSL Gnu Scientific Library          https://www.gnu.org/software/gsl/
include(cmake/Find_or_install_SPECTRA.cmake)             ### For Spectra eigenvalue solver           https://github.com/yixuan/spectra





message("======================================")
message("=== SUMMARY OF CMAKE CONFIGURATION ===")
message("======================================")
message("FOUND EIGEN3   :   ${EIGEN3_INCLUDE_DIR}")
message("FOUND HDF5     :   ${HDF5_LIBRARIES}")
message("FOUND BLAS     :   ${BLAS_LIBRARIES}")
message("FOUND LAPACK   :   ${LAPACK_LIBRARIES}")
message("FOUND ARPACK   :   ${ARPACK_LIBRARIES}")
message("FOUND ARPACKPP :   ${ARPACKPP_INCLUDE_DIR}")
message("FOUND GSL      :   ${GSL_LIBRARIES}")
message("FOUND SPECTRA  :   ${SPECTRA_INCLUDE_DIR}")
message("FOUND INTEL MKL:   ${MKL_LIBRARIES}")
message("======================================")


###################
### INFORMATION ###
###################
message("==================================================================")
message("")
message("      To clean cmake files, run './build.sh clean'")
message("      To clean cmake files downloaded libraries, run './clean.sh '")
message("      To launch this program,  './run.sh'")
message("")
message("==================================================================")





# To print all variables, use the code below:
#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
