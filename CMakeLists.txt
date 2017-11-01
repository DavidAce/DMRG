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



set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if (NOT "${GSL_DIR}" MATCHES "hunter")
    find_package(GSL)                                       ### Find and define includes for Eigen Library
    message("Searching for native GSL installation. Found: ${GSL_FOUND}")
endif()
if (NOT GSL_FOUND)
    message("Using GSL provided by Hunter")
        hunter_add_package(GSL)                             ### Download and install GSL Library
        find_package(GSL CONFIG REQUIRED)                   ### Find and define includes for GSL Library
endif()



if (NOT "${hdf5_DIR}" MATCHES "hunter")
    find_package(HDF5 1.8 MODULE COMPONENTS CXX)   ### Find HDF5 installed in system
#    include_directories(${HDF5_INCLUDE_DIRS})
    add_definitions    (${HDF5_DEFINITIONS})
#    link_libraries     (${HDF5_LIBRARIES})
    message("Searching for native HDF5 installation. Found: ${HDF5_FOUND}")
endif()
if(NOT HDF5_FOUND)
    message("Using HDF5 provided by Hunter")
    hunter_add_package(hdf5)                            ### Download and install HDF5 Library
    find_package(ZLIB CONFIG REQUIRED)                  ### Require ZLIB (for HDF5)
    find_package(szip CONFIG REQUIRED)                  ### Require szip (for HDF5)
    find_package(hdf5 CONFIG REQUIRED COMPONENTS CXX)   ### Require HDF5
endif()


if (NOT "${Eigen3_DIR}" MATCHES "hunter")
    find_package(Eigen3 3.3)                            ### Find and define includes for Eigen Library
    message("Searching for native Eigen3 installation. Found: ${HDF5_FOUND}")
endif()

if (NOT EIGEN3_FOUND)
    message("Using Eigen3 provided by Hunter")
    hunter_add_package(Eigen)                           ### Download and install Eigen library
    find_package(Eigen3 CONFIG REQUIRED)                ### Find and define includes for Eigen Library
endif()


##
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

###################
### Setup flags ###
###################
set(COMMON_OPTIONS -g -Wall -march=native)
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    message("Special flags for g++:  -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context")
    set(COMMON_OPTIONS ${COMMON_OPTIONS}   -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Special flags for Clang: -Wno-invalid-partial-specialization -Wno-missing-braces")
    set(COMMON_OPTIONS ${COMMON_OPTIONS} -Wno-invalid-partial-specialization -Wno-missing-braces)
endif()


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


target_include_directories(DMRG++ PRIVATE source)                   ### Project source filess
target_include_directories(DMRG++ PRIVATE spectra/include)          ### Spectra                                 https://spectralib.org/
target_include_directories(DMRG++ PRIVATE ${HDF5_INCLUDE_DIRS})     ### Spectra                                 https://spectralib.org/
target_link_libraries(DMRG++ -lstdc++fs)                            ### For experimental::filesystem
target_link_libraries(DMRG++ GSL::gsl)                              ### For GSL GNU Scientific Library          https://www.gnu.org/software/gsl/
target_link_libraries(DMRG++ ${HDF5_LIBRARIES})                     ### For HDF5                                https://support.hdfgroup.org/HDF5/
if("${Eigen3_DIR}" MATCHES "hunter")
    target_link_libraries(DMRG++ Eigen3::Eigen)                     ### For Eigen linear algebra header-library http://eigen.tuxfamily.org/
else()
    target_include_directories(DMRG++ PRIVATE ${EIGEN3_INCLUDE_DIR})### For Eigen linear algebra header-library http://eigen.tuxfamily.org/
endif()