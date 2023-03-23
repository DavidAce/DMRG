if(DMRG_PACKAGE_MANAGER MATCHES "find")
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    #    find_package(Fortran REQUIRED)
    ##############################################################################
    ###  Optional Intel MKL support. Uses OpenBLAS as fall-back                ###
    ##############################################################################
    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS blas lapack gf gnu_thread lp64 REQUIRED MODULE)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    else()
        find_package(OpenBLAS 0.3.20 REQUIRED)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
        #For convenience, define these targes
        add_library(BLAS::BLAS ALIAS OpenBLAS::OpenBLAS)
        add_library(LAPACK::LAPACK ALIAS OpenBLAS::OpenBLAS)
        add_library(lapacke::lapacke ALIAS OpenBLAS::OpenBLAS)
    endif()


    find_package(Lapacke           REQUIRED)                                         # Lapacke needed by arpack++, included in MKL or OpenBLAS
    find_package(pcg-cpp           REQUIRED)
    find_package(Eigen3     3.4.0  REQUIRED)                                         # Eigen3 numerical library (needed by ceres and h5pp)
    find_package(h5pp       1.11.0 REQUIRED)                                         # h5pp for writing to file binary in format
    find_package(Ceres      2.0    REQUIRED PATH_SUFFIXES ceres ceres/lib)           # ceres-solver (for L-BFGS routine)
    find_package(CLI11      2.1.1  REQUIRED)                                         # Command line argument parser
    find_package(arpack-ng  3.8.0  REQUIRED)                                         # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    find_package(libunwind  1.6.2  REQUIRED)
    find_package(Backward   1.6    REQUIRED)
    # Arpack needs to link to extra libraries before moving on
    if(arpack-ng_FOUND AND TARGET ARPACK::ARPACK)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()
    if(TARGET arpack-ng::arpack-ng)
        add_library(ARPACK::ARPACK ALIAS arpack-ng::arpack-ng)
    endif()
    find_package(arpack++ 2.3.0 REQUIRED)                                                 # C++ frontend for arpack-ng. Custom find module.


    target_link_libraries(dmrg-deps INTERFACE
            CLI11::CLI11
            pcg-cpp::pcg-cpp
            h5pp::h5pp
            arpack++::arpack++
            primme::primme
            Ceres::ceres
            Backward::Backward
            libunwind::libunwind
            BLAS::BLAS

            )

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
endif()