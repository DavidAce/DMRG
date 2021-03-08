if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")

    ##############################################################################
    ###  Optional OpenMP support                                               ###
    ###  Note that Clang has some  trouble with static openmp and that         ###
    ###  and that static openmp is not recommended. This tries to enable       ###
    ###  static openmp anyway because I find it useful. Installing             ###
    ###  libiomp5 might help for shared linking.                               ###
    ##############################################################################
    if(DMRG_ENABLE_OPENMP)
        find_package(OpenMP) # Uses DMRG's own find module
    endif()
    find_package(Fortran REQUIRED)
    include(cmake-modules/SetupMKL.cmake)                           # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    include(cmake-modules/Get_OpenBLAS.cmake)                     # If MKL is not on openblas will be used instead. Includes lapack.
    include(cmake-modules/Get_Eigen3.cmake)                       # Eigen3 numerical library (needed by ceres and h5pp)
    include(cmake-modules/Get_h5pp.cmake)                         # h5pp for writing to file binary in format
    include(cmake-modules/Get_arpack-ng.cmake)                    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    include(cmake-modules/Get_arpack++.cmake)                     # C++ frontend for arpack-ng
    include(cmake-modules/Get_gflags.cmake)                       # Google Flags library needed by ceres-solver
    include(cmake-modules/Get_glog.cmake)                         # Google logging library needed by ceres-solver
    include(cmake-modules/Get_ceres-solver.cmake)                 # ceres-solver (for L-BFGS routine)


    if(TARGET Eigen3::Eigen AND DMRG_ENABLE_THREADS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    endif()

    if(TARGET Eigen3::Eigen AND TARGET BLAS::BLAS )
        if(TARGET mkl::mkl)
            message(STATUS "Eigen3 will use MKL")
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
        else ()
            message(STATUS "Eigen3 will use OpenBLAS")
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE BLAS::BLAS)
        endif()

        # AVX2 aligns 32 bytes (AVX512 aligns 64 bytes).
        # When running on Tetralith, with march=native, there can be alignment mismatch
        # in ceres which results in a segfault on free memory.
        # Something like "double free or corruption ..."
        #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith
        cmake_host_system_information(RESULT _host_name  QUERY HOSTNAME)
        if(_host_name MATCHES "tetralith|triolith")
            message(STATUS "Applying special Eigen compile definitions for Tetralith: EIGEN_MAX_ALIGN_BYTES=16")
            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=1) # May work to fix CERES segfaults!!!
            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=32)
        else()
            message(STATUS "Applying special Eigen compile definitions for general machines")
#            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=1) # May work to fix CERES segfaults!!!
#            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=32)
        endif()

    endif()


    ##################################################################
    ### Link all the things!                                       ###
    ##################################################################
    if(TARGET h5pp::h5pp)
        target_link_libraries(dmrg-main PUBLIC h5pp::h5pp)
    endif()
    if(TARGET spdlog::spdlog)
        target_link_libraries(dmrg-opt PUBLIC spdlog::spdlog)
        target_link_libraries(dmrg-eig PUBLIC spdlog::spdlog)
        target_link_libraries(dmrg-arp PUBLIC spdlog::spdlog)
    endif()
    if(TARGET Eigen3::Eigen)
        target_link_libraries(dmrg-main PUBLIC Eigen3::Eigen)
        target_link_libraries(dmrg-eig PUBLIC Eigen3::Eigen)
        target_link_libraries(dmrg-opt PUBLIC Eigen3::Eigen)
    endif()
    if(TARGET Ceres::ceres)
        target_link_libraries(dmrg-opt PUBLIC Ceres::ceres)
    endif()
    if(TARGET ARPACK::ARPACK++)
        target_link_libraries(dmrg-arp PUBLIC ARPACK::ARPACK++)
    endif()
    if(TARGET openmp::openmp)
        target_link_libraries(dmrg-flags INTERFACE openmp::openmp)
    else()
        target_compile_options(dmrg-flags INTERFACE -Wno-unknown-pragmas)
    endif()
    if(TARGET unwind::unwind)
        target_link_libraries(dmrg-dbg PUBLIC unwind::unwind)
        target_compile_definitions(dmrg-dbg PUBLIC DMRG_HAS_UNWIND=1)
    endif()
    set(DMRG_TARGETS h5pp::h5pp ARPACK::ARPACK++ Ceres::ceres Eigen3::Eigen dmrg-flags)
    mark_as_advanced(DMRG_TARGETS)
endif()