cmake_minimum_required(VERSION 3.7)

### Hunter package manager ###
### For newest releases go to
###     https://github.com/ruslo/hunter/releases
include("cmake/HunterGate.cmake")
HunterGate(
        URL "https://github.com/ruslo/hunter/archive/v0.19.157.tar.gz"
        SHA1 "4bd7cbd1fb5dd33f52e19e05d893d362712cea9e"
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
    set(GSL_FOUND 1)
    message("Using GSL provided by Hunter")
    hunter_add_package(GSL)                             ### Download and install GSL Library
    find_package(GSL CONFIG REQUIRED)                   ### Find and define includes for GSL Library
endif()



if (NOT "${hdf5_DIR}" MATCHES "hunter")
    find_package(HDF5 1.8 MODULE COMPONENTS CXX)   ### Find HDF5 installed in system
    add_definitions    (${HDF5_DEFINITIONS})
    message("Searching for native HDF5 installation. Found: ${HDF5_FOUND}")
endif()
if(NOT HDF5_FOUND)
    set(HDF5_FOUND 1)
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
    set(EIGEN3_FOUND 1)
    message("Using Eigen3 provided by Hunter")
    hunter_add_package(Eigen)                           ### Download and install Eigen library
    find_package(Eigen3 CONFIG REQUIRED)                ### Find and define includes for Eigen Library
endif()




###################
### Setup flags ###
###################
set(COMMON_OPTIONS -g -Wall -march=native )
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3 -flto)
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
target_link_libraries(DMRG++ -flto)                                             ### Link time optimization flag

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

#
#if (NOT "${BLAS_DIR}" MATCHES "hunter")
#    find_package(BLAS)                                      ### Find and define includes for Eigen Library
#    message("Searching for native OpenBLAS installation. Found: ${BLAS_FOUND}")
#endif()
#if (NOT BLAS_FOUND)
#    set(BLAS_FOUND 1)
#    message("Using OpenBLAS provided by Hunter")
#    hunter_add_package(OpenBLAS)                             ### Download and install GSL Library
#    find_package(OpenBLAS CONFIG REQUIRED)                   ### Find and define includes for GSL Library
#endif()

#
#include(cmake/DownloadProject.cmake)
#download_project(PROJ                arpackng
#        GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
#        GIT_TAG             master
#        SOURCE_DIR "${PROJECT_BINARY_DIR}/deps/arpack"
#        BINARY_DIR "${PROJECT_BINARY_DIR}/deps/arpack-build"
#        INSTALL_DIR "${PROJECT_BINARY_DIR}/deps/arpack-install"
#        GIT_PROGRESS 1
#        CMAKE_ARGS
#        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#        -DINSTALL_BIN_DIR:PATH=<INSTALL_DIR>/bin
#        -DINSTALL_INC_DIR:PATH=<INSTALL_DIR>/include
#        -DINSTALL_LIB_DIR:PATH=<INSTALL_DIR>/lib
#        -DINSTALL_MAN_DIR:PATH=<INSTALL_DIR>/share/man
#        -DINSTALL_PKGCONFIG_DIR:PATH=<INSTALL_DIR>/share/pkgconfig
#        )
#
#add_subdirectory(${arpackng_SOURCE_DIR} ${arpackng_BINARY_DIR})
#


#if("${OpenBLAS_DIR}" MATCHES "hunter")
#target_link_libraries(DMRG++ OpenBLAS::OpenBLAS)
#endif()
#target_link_libraries(DMRG++ arpackng)


#include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
#
#if (NOT ARPACK_FOUND)
#    message(STATUS "ARPACK has not been found. Therefore it will automatically be downloaded during the compilation process.")
#    ExternalProject_Add(ARPACK
#            SOURCE_DIR "${PROJECT_BINARY_DIR}/deps/arpack"
#            BINARY_DIR "${PROJECT_BINARY_DIR}/deps/arpack-build"
#            INSTALL_DIR "${PROJECT_BINARY_DIR}/deps/arpack-install"
#            GIT_REPOSITORY "${git_protocol}://github.com/opencollab/arpack-ng.git"
#            GIT_PROGRESS 1
#            CMAKE_ARGS
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DINSTALL_BIN_DIR:PATH=<INSTALL_DIR>/bin
#            -DINSTALL_INC_DIR:PATH=<INSTALL_DIR>/include
#            -DINSTALL_LIB_DIR:PATH=<INSTALL_DIR>/lib
#            -DINSTALL_MAN_DIR:PATH=<INSTALL_DIR>/share/man
#            -DINSTALL_PKGCONFIG_DIR:PATH=<INSTALL_DIR>/share/pkgconfig)
#    ExternalProject_Get_Property(ARPACK install_dir)
#    set(ARPACK_LIB_DIR ${install_dir}/lib)
#    set(ARPACK_INC_DIR ${install_dir}/include)
#    set(ARPACK_BIN_DIR ${install_dir}/bin)
#    add_library(arpack STATIC IMPORTED)
#    set_target_properties(ARPACK PROPERTIES IMPORTED_LOCATION ${ARPACK_LIB_DIR}/arpack.lib)
#
#endif()


#set(${MKLROOT} /opt)
#target_include_directories(DMRG++ PUBLIC ${MKL_ROOT}/include)

#target_link_libraries(DMRG++  )


#message("${LIBRARY_MKL}")
#message(STATUS "LIBRARY_MKL = " $ENV{LIBRARY_MKL})
#
#if("${MKLROOT}" MATCHES "intel")
#    message("Found Intel MKL")
#    set(COMMON_OPTIONS ${COMMON_OPTIONS}  -m64 -I${MKLROOT}/include )
#endif()
#set(MKLROOT /opt/intel/compilers_and_libraries_2018/linux/mkl)
#message("MKLROOT = " $ENV{MKLROOT})
#message("MKLROOT = " ${MKLROOT})
#set(COMMON_OPTIONS ${COMMON_OPTIONS} -DMKL_LP64 -m64 -I${MKLROOT}/include  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl)





# To print all variables, use the code below:
#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()