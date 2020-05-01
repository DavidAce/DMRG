if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    #  Make sure we use DMRG's own find modules
    list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake-modules)
    # Make sure find_package looks for the packages we just built instead of rebuilding them
    # when we rerun the cmake config.
    list(INSERT CMAKE_PREFIX_PATH 0 ${CMAKE_INSTALL_PREFIX})

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
    include(cmake-modules/Find_dont_install_INTELMKL.cmake)         # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    include(cmake-modules/Fetch_OpenBLAS.cmake)                     # If MKL is not on openblas will be used instead. Includes lapack.
    include(cmake-modules/Fetch_arpack-ng.cmake)                    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    include(cmake-modules/Fetch_arpack++.cmake)                     # C++ frontend for arpack-ng
    include(cmake-modules/Fetch_Eigen3.cmake)                       # Eigen3 numerical library (needed by ceres and h5pp)
    include(cmake-modules/Fetch_gflags.cmake)                       # Google Flags library needed by ceres-solver
    include(cmake-modules/Fetch_glog.cmake)                         # Google logging library needed by ceres-solver
    include(cmake-modules/Fetch_ceres-solver.cmake)                 # ceres-solver (for L-BFGS routine)
    include(cmake-modules/Fetch_h5pp.cmake)                         # h5pp for writing to file binary in format


    ##################################################################
    ### Link all the things!                                       ###
    ##################################################################
    if(TARGET ceres::ceres)
        list(APPEND NATIVE_TARGETS ceres::ceres)
    endif()
    if(TARGET Eigen3::Eigen)
        list(APPEND NATIVE_TARGETS Eigen3::Eigen)
    endif()
    if(TARGET h5pp::h5pp)
        list(APPEND NATIVE_TARGETS h5pp::h5pp)
    endif()
    if(TARGET arpack::arpack++)
        list(APPEND NATIVE_TARGETS arpack::arpack++)
    endif()
    if(TARGET openmp::openmp)
        list(APPEND NATIVE_TARGETS openmp::openmp)
    else()
        target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
    endif()
    if(TARGET Threads::Threads)
        list(APPEND NATIVE_TARGETS Threads::Threads)
    endif()
    if(NATIVE_TARGETS)
        mark_as_advanced(NATIVE_TARGETS)
    endif()



    if(TARGET Eigen3::Eigen AND TARGET openmp::openmp)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_THREADS)
    endif()

    if(TARGET Eigen3::Eigen AND TARGET blas::blas )
        set(EIGEN3_USING_BLAS ON)
        if(TARGET mkl::mkl)
            message(STATUS "Eigen3 will use MKL")
            target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_MKL_ALL)
            target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
        else ()
            message(STATUS "Eigen3 will use OpenBLAS")
            target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_BLAS)
            target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE blas::blas)
        endif()
        # Use this flag if Ceres is giving you trouble!
        # For some reason it starts mixing aligned and hand-made aligned malloc and freeing them willy nilly
        # This flag forces its hand and avoids a segfault in some cases.
        #    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_MALLOC_ALREADY_ALIGNED=0) # Finally something works!!!
    endif()
endif()