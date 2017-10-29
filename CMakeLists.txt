cmake_minimum_required(VERSION 3.7)

### Hunter package manager ###
### For newest releases go to
###     https://github.com/ruslo/hunter/releases
include("cmake/HunterGate.cmake")
HunterGate(
        URL "https://github.com/ruslo/hunter/archive/v0.19.139.tar.gz"
        SHA1 "31a0f9f2fb1b7dc44c7c19656e2a8acabe284959"
)


project(DMRG++)
set(CMAKE_VERBOSE_MAKEFILE OFF)

################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake/gitversion.cmake)


hunter_add_package(Eigen)               ### Download and install Eigen library
find_package(Eigen3 CONFIG REQUIRED)    ### Find and define includes for Eigen Library

hunter_add_package(GSL)                 ### Download and install GSL Library
find_package(GSL CONFIG REQUIRED)       ### Find and define includes for GSL Library

hunter_add_package(hdf5)                ### Download and install HDF5 Library
find_package(ZLIB CONFIG REQUIRED)      ### Require ZLIB (for HDF5)
find_package(szip CONFIG REQUIRED)      ### Require szip (for HDF5)
find_package(hdf5 CONFIG REQUIRED)      ### Require HDF5
find_package(OpenMP REQUIRED)           ### Require OpenMP


###################
### Setup flags ###
###################
set(COMMON_OPTIONS -g -fopenmp -Wall -march=native -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context)
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(COMMON_OPTIONS ${COMMON_OPTIONS}  -Wno-invalid-partial-specialization)
endif()
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")


###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(DMRG++ ${SOURCES} ${HEADERS})


########################################################
###  target properties go AFTER add_executable()   ###
########################################################
target_compile_features(DMRG++ PUBLIC cxx_std_17)                               ### Demand c++17 compiler support
target_compile_options(DMRG++ PUBLIC ${COMMON_OPTIONS})                         ### Common options
target_compile_options(DMRG++ PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")     ### Debug build options
target_compile_options(DMRG++ PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>") ### Release build options


target_include_directories(DMRG++ PRIVATE source)           ### Project source filess
target_include_directories(DMRG++ PRIVATE spectra/include)  ### Spectra                                 https://spectralib.org/
target_link_libraries(DMRG++ -lstdc++fs -fopenmp)           ### For experimental::filesystem and OpenMP
target_link_libraries(DMRG++ Eigen3::Eigen)                 ### For Eigen linear algebra header-library http://eigen.tuxfamily.org/
target_link_libraries(DMRG++ GSL::gsl)                      ### For GSL GNU Scientific Library          https://www.gnu.org/software/gsl/
target_link_libraries(DMRG++ hdf5)                          ### For HDF5                                https://support.hdfgroup.org/HDF5/

