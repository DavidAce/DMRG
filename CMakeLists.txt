cmake_minimum_required(VERSION 3.8)
enable_language(Fortran)
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



###########################################
###  Find or install all dependencies   ###
###########################################
include(cmake/Find_or_install_EIGEN3.cmake)                                     ### For Eigen linear algebra header-library http://eigen.tuxfamily.org/
include(cmake/Find_or_install_HDF5.cmake)                                       ### For HDF5 data env_storage format            https://support.hdfgroup.org
include(cmake/Find_or_install_LAPACK_and_BLAS.cmake)                            ### For Lapack and Blas libraries           https://github.com/Reference-LAPACK/lapack
include(cmake/Find_or_install_ARPACK.cmake)                                     ### For Arpack-ng, a mantained arpack fork  https://github.com/opencollab/arpack-ng
include(cmake/Find_or_install_ARPACKPP.cmake)                                   ### For ARPACKPP header-arplibrary          https://github.com/m-reuter/arpackpp
include(cmake/Find_or_install_GSL.cmake)                                        ### For GSL Gnu Scientific Library          https://www.gnu.org/software/gsl/
include(cmake/Find_or_install_SPECTRA.cmake)                                    ### For Spectra eigenvalue solver           https://github.com/yixuan/spectra
include(cmake/FindGFortran.cmake)                                               ### For Fortran library

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
message("======================================")

include_directories(libs/spectra/include)



###################
### Setup flags ###
###################
set(COMMON_OPTIONS -g -Wall -march=native)
set(DEBUG_OPTIONS       -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    message("Special flags for g++:  -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-maybe-uninitialized")
    set(COMMON_OPTIONS ${COMMON_OPTIONS}  -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-maybe-uninitialized)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Special flags for Clang:  -stdlib=libstdc++ -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-invalid-partial-specialization -Wno-missing-braces -Wno-overloaded-virtual -Wno-maybe-uninitialized")
    set(COMMON_OPTIONS ${COMMON_OPTIONS}   -stdlib=libstdc++ -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-invalid-partial-specialization -Wno-missing-braces -Wno-overloaded-virtual -Wno-maybe-uninitialized)
endif()



###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(DMRG++ ${SOURCES} ${HEADERS})


###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_features(DMRG++ PUBLIC cxx_std_17)                                         ### Demand c++17 compiler support
target_compile_options(DMRG++ PUBLIC ${COMMON_OPTIONS})                                   ### Common options
target_compile_options(DMRG++ PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options(DMRG++ PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_link_libraries(DMRG++ -flto)                                                       ### Link time optimization flag


#######################################
###  Link and include dependencies  ###
#######################################

target_include_directories(DMRG++ PRIVATE source)
target_include_directories(DMRG++ PRIVATE ${EIGEN3_INCLUDE_DIR})
target_include_directories(DMRG++ PRIVATE ${ARPACKPP_INCLUDE_DIR})
target_include_directories(DMRG++ PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(DMRG++ PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(DMRG++ PRIVATE libs/spectra/include)
target_link_libraries(DMRG++ -lstdc++fs)                              ### For experimental::filesystem
target_link_libraries(DMRG++ ${ARPACK_LIBRARIES})
target_link_libraries(DMRG++ ${LAPACK_LIBRARIES})
target_link_libraries(DMRG++ ${BLAS_LIBRARIES})
target_link_libraries(DMRG++ ${HDF5_LIBRARIES} -ldl)
target_link_libraries(DMRG++ ${GSL_LIBRARIES})
target_link_libraries(DMRG++ ${GFORTRAN_LIB})







# To print all variables, use the code below:
#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
