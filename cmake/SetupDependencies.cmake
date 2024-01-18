
##############################################################################
###  Optional Intel MKL support. Uses OpenBLAS as fall-back                ###
##############################################################################
if(NOT DEFINED BLA_VENDOR AND NOT DEFINED $ENV{BLA_VENDOR})
    message(FATAL_ERROR  "BLA_VENDOR is undefined: https://cmake.org/cmake/help/latest/module/FindBLAS.html" )
endif()
if(BLA_VENDOR MATCHES Intel OR ENV{BLA_VENDOR} MATCHES Intel)
    find_package(MKL REQUIRED BYPASS_PROVIDER) # Will try to find MKLconfig.cmake in the Intel MKL install dir
endif()



# Setup dependencies
find_package(Threads                    REQUIRED)
find_package(OpenMP                     REQUIRED COMPONENTS CXX )
find_package(gfortran                   REQUIRED OPTIONAL_COMPONENTS quadmath)
find_package(Lapacke                    REQUIRED MODULE)
find_package(pcg-cpp                    REQUIRED)
find_package(Eigen3     3.4.0           REQUIRED)                                         # Eigen3 numerical library (needed by ceres and h5pp)
find_package(h5pp       1.11.0...1.11.1 REQUIRED)                                         # h5pp for writing to file binary in format
find_package(fmt        10.1.0...10.1.2 REQUIRED)
find_package(spdlog     1.11.0...1.12.0 REQUIRED)
find_package(CLI11      2.1.1...2.3.2   REQUIRED)                                         # Command line argument parser
find_package(Backward   1.6             REQUIRED)
#find_package(arpack++   2.3.0  REQUIRED)                                          # C++ frontend for arpack-ng. Custom find module.
#find_package(mpfr       4.1.0  REQUIRED)

include(cmake/CheckCompile.cmake)
check_compile(Lapacke lapacke::lapacke REQUIRED)


# Install dependencies that need manual installation
include(cmake/cmake_dependency_provider/PKGInstall.cmake)
pkg_install(arpack-ng)
pkg_install(arpack++)
pkg_install(primme)

find_package(arpack-ng 3.8.0...3.9.0 REQUIRED MODULE BYPASS_PROVIDER)
find_package(arpack++                REQUIRED MODULE BYPASS_PROVIDER)
find_package(primme                  REQUIRED MODULE BYPASS_PROVIDER)

# Link all dependencies to dmrg-deps
if(NOT TARGET dmrg-deps)
    add_library(dmrg-deps INTERFACE)
endif()
if(NOT TARGET dmrg-flags)
    add_library(dmrg-flags INTERFACE)
endif()
target_link_libraries(dmrg-flags INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(dmrg-deps INTERFACE
            CLI11::CLI11
            pcg-cpp::pcg-cpp
            h5pp::h5pp
            Eigen3::Eigen
            fmt::fmt
            spdlog::spdlog
            lapacke::lapacke
            arpack++::arpack++
            arpack-ng::arpack-ng
            primme::primme
            # We link Backward::Backward on the dmrg-stacktrace object directly
            )

if(DMRG_ENABLE_TBLIS)
    pkg_install(tblis)
    find_package(tblis REQUIRED MODULE BYPASS_PROVIDER)
    target_link_libraries(dmrg-deps INTERFACE tblis::tblis)
endif()

# Configure Eigen
if(TARGET Eigen3::Eigen)
    if(EIGEN_USE_THREADS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    endif()
    if(TARGET MKL::MKL)
        message(STATUS "Eigen will use MKL")
        target_link_libraries(Eigen3::Eigen INTERFACE MKL::MKL)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
    elseif(TARGET BLAS::BLAS)
        message(STATUS "Eigen will use BLAS")
        target_link_libraries(Eigen3::Eigen INTERFACE BLAS::BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE)
    endif()
else()
    message(FATAL_ERROR "Target not defined: Eigen3::Eigen")
endif()

### Set the floating point type high-precision arithmetic (used in lbit Hamiltonian parameters
### for accurate long time-scale evolution)
if(DMRG_USE_QUADMATH)
    find_package(quadmath REQUIRED)
    target_compile_definitions(dmrg-flags INTERFACE USE_QUADMATH)
    target_link_libraries(dmrg-flags INTERFACE quadmath::quadmath)
endif ()