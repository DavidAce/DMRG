
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
find_package(Threads           REQUIRED)
find_package(OpenMP            REQUIRED COMPONENTS CXX )
find_package(gfortran          REQUIRED COMPONENTS quadmath)
find_package(Lapacke           REQUIRED MODULE)
find_package(pcg-cpp           REQUIRED)
find_package(Eigen3     3.4.0  REQUIRED)                                         # Eigen3 numerical library (needed by ceres and h5pp)
find_package(h5pp       1.11.1 REQUIRED)                                         # h5pp for writing to file binary in format
find_package(fmt        9.1.0  REQUIRED)
find_package(spdlog     1.10.0 REQUIRED)
find_package(Ceres      2.0    REQUIRED)                                         # ceres-solver (for L-BFGS routine)
find_package(CLI11      2.1.1  REQUIRED)                                         # Command line argument parser
find_package(arpack-ng  3.8.0  REQUIRED)                                         # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
find_package(Backward   1.6    REQUIRED)
find_package(arpack++   2.3.0  REQUIRED)                                          # C++ frontend for arpack-ng. Custom find module.
#find_package(mpfr       4.1.0  REQUIRED)

include(cmake/CheckCompile.cmake)
check_compile(Lapacke lapacke::lapacke REQUIRED)


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
            arpack++::arpack++
            primme::primme
            Ceres::ceres
            lapacke::lapacke
#            mpfr::mpfr
            # We link Backward::Backward on the dmrg-stacktrace object directly
            )


# Install dependencies that need manual installation
include(cmake/cmake_dependency_provider/PKGInstall.cmake)


#pkg_install(mpreal) # For MPFRC++ - c++ frontend for the mpfr multiprecision library
#find_package(mpreal REQUIRED MODULE)
pkg_install(primme)
find_package(primme REQUIRED MODULE)
target_link_libraries(dmrg-deps INTERFACE primme::primme)

if(DMRG_ENABLE_TBLIS)
    pkg_install(tblis)
    find_package(tblis REQUIRED MODULE)
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

# Fix issue with Ceres linking to cuda
find_package(CUDA) # Same call as when building Ceres
if (CUDA_FOUND)
    message("-- Found CUDA version ${CUDA_VERSION}: "
            "${CUDA_LIBRARIES};"
            "${CUDA_cusolver_LIBRARY};"
            "${CUDA_cusparse_LIBRARY};"
            "${CUDA_CUBLAS_LIBRARIES}"
            )
    target_link_libraries(dmrg-deps INTERFACE ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
else ()
    target_compile_definitions(dmrg-deps INTERFACE CERES_NO_CUDA)
endif ()


### Set the floating point type high-precision arithmetic (used in lbit Hamiltonian parameters
### for accurate long time-scale evolution)
if(DMRG_USE_QUADMATH)
    find_package(quadmath REQUIRED)
    target_compile_definitions(dmrg-flags INTERFACE USE_QUADMATH)
    target_link_libraries(dmrg-flags INTERFACE quadmath::quadmath)
endif ()