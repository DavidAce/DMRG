if(DMRG_PACKAGE_MANAGER MATCHES "find")
    if(NOT BUILD_SHARED_LIBS)
        set(GFLAGS_COMPONENTS COMPONENTS)
        set(GFLAS_ITEMS nothreads_static)
    endif()
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    find_package(Fortran REQUIRED)

    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS gf gnu_thread lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    else()
        find_package(OpenBLAS 0.3.8 REQUIRED) # If MKL is not on openblas will be used instead. Includes lapack.
        if(TARGET OpenBLAS::OpenBLAS)
            target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
            #For convenience, define these targes
            add_library(BLAS::BLAS ALIAS OpenBLAS::OpenBLAS)
            add_library(LAPACK::LAPACK ALIAS OpenBLAS::OpenBLAS)
            add_library(lapacke::lapacke  ALIAS OpenBLAS::OpenBLAS)
        endif()

    endif()
    find_package(Lapacke          REQUIRED)                                         # Lapacke needed by arpack++, included in MKL or OpenBLAS
    find_package(Eigen3     3.3.7 REQUIRED)                                         # Eigen3 numerical library (needed by ceres and h5pp)
    find_package(h5pp       1.9.0 REQUIRED)                                         # h5pp for writing to file binary in format
    find_package(gflags     2.2.2 REQUIRED ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS})    # Google Flags library needed by ceres-solver
    find_package(glog       0.4   REQUIRED)                                         # Google logging library needed by ceres-solver
    find_package(Ceres      2.0   REQUIRED PATH_SUFFIXES ceres ceres/lib)           # ceres-solver (for L-BFGS routine)
    find_package(CLI11      2.1.1 REQUIRED)                                         # Command line argument parser
    find_package(arpack-ng  3.8.0 REQUIRED)                                         # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    # Arpack needs to link to extra libraries before moving on
    if(arpack-ng_FOUND)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()
    find_package(arpack++ REQUIRED)                                                 # C++ frontend for arpack-ng. Custom find module.
endif()