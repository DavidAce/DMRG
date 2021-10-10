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
    endif()
    find_package(Lapacke REQUIRED)                                   # Lapacke needed by arpack++, included in MKL or OpenBLAS

    find_package(Eigen3 3.3.7 REQUIRED)      # Eigen3 numerical library (needed by ceres and h5pp)
    find_package(h5pp 1.9.0 REQUIRED)        # h5pp for writing to file binary in format
    find_package(arpack-ng 3.8.0 REQUIRED)   # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.

    # Arpack needs to link to extra libraries before moving on
    if(arpack-ng_FOUND)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()

    find_package(arpack++ REQUIRED)                                  # C++ frontend for arpack-ng. Custom find module.
    if(arpack++_FOUND)
        target_link_libraries(arpack::arpack++ INTERFACE ARPACK::ARPACK lapacke::lapacke)
    endif()
    find_package(gflags 2.2.2 ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS} REQUIRED)      # Google Flags library needed by ceres-solver
    find_package(glog 0.4 REQUIRED)          # Google logging library needed by ceres-solver
    find_package(Ceres 2.0 PATH_SUFFIXES ceres ceres/lib REQUIRED)      # ceres-solver (for L-BFGS routine)
    find_package(CLI11 2.1.1 REQUIRED)

    if(TARGET Eigen3::Eigen AND DMRG_ENABLE_THREADS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    endif()

    if(TARGET Eigen3::Eigen)
        if(MKL_FOUND)
            message(STATUS "Eigen3 will use MKL")
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
        elseif(OpenBLAS_FOUND)
            message(STATUS "Eigen3 will use OpenBLAS")
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE OpenBLAS::OpenBLAS)
        endif()

        cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
        if(_host_name MATCHES "tetralith|triolith")
            # AVX aligns 32 bytes (AVX512 aligns 64 bytes).
            # When running on Tetralith, with march=native, there can be alignment mismatch
            # in ceres which results in a segfault on free memory.
            # Something like "double free or corruption ..."
            #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith

            ### NOTE October 4 2020 ####
            #
            # Another flag that seems to fix weird release-only bugs is
            #       -fno-strict-aliasing

            ### NOTE August 26 2020 ####
            #
            # Ceres started crashing on Tetralith again using -march=native.
            # Tried to solve this issue once and for all.
            # I've tried the following flags during compilation of DMRG++ and ceres-solver:
            #
            #           -DEIGEN_MALLOC_ALREADY_ALIGNED=[none,0,1]
            #           -DEIGEN_MAX_ALIGN_BYTES=[none,16,32]
            #           -march=[none,native]
            #           -std=[none,c++17]
            #
            # Up until now, [0,16,none,none] has worked but now for some reason it stopped now.
            # I noticed the stdc++=17 flag was not being passed on conan builds, so ceres defaulted to -std=c++14 instead.
            # I fixed this in the conanfile.py of the ceres build. The -package-manager=cmake method already had this fixed.
            # When no Eigen flags were passed, and ceres-solver finally built with -std=c++17 the issues vanished.
            # In the end what worked was [none,none,native,c++17] in both DMRG++ and ceres-solver.
            # It is important that the same eigen setup is used in all compilation units, and c++17/c++14 seems to
            # make Eigen infer some of the flags differently. In any case, settinc c++17 and no flags for eigen anywhere
            # lets Eigen do its thing in the same way everywhere.

            #            message(STATUS "Applying special Eigen compile definitions for Tetralith: EIGEN_MAX_ALIGN_BYTES=16")
            #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfault?
            #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
        else()
            #            message(STATUS "Applying special Eigen compile definitions for general machines: EIGEN_MAX_ALIGN_BYTES=16")
            #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfaults!!!
            #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
        endif()
    else()
        message(FATAL_ERROR "Target not defined: Eigen3::Eigen")
    endif()

endif()