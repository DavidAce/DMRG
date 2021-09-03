
if (DMRG_PACKAGE_MANAGER STREQUAL "cmake")
    include(cmake/InstallPackage.cmake)

    list(APPEND h5pp_ARGS -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND h5pp_ARGS -DH5PP_PACKAGE_MANAGER:STRING=cmake)
    list(APPEND h5pp_ARGS -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})


    list(APPEND Ceres_ARGS -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND Ceres_ARGS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND Ceres_ARGS -Dglog_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})

    if (NOT BUILD_SHARED_LIBS)
        set(GFLAGS_COMPONENTS COMPONENTS)
        set(GFLAS_ITEMS nothreads_static)
    endif ()

    # Find packages or install if missing
    find_package(Threads REQUIRED)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    find_package(Fortran REQUIRED)

    if (DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS blas lapack gf gnu_thread lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    endif ()

    if (NOT MKL_FOUND)
        # If MKL is not on openblas will be used instead. Includes blas, lapack and lapacke
        if(NOT DEFINED OPENBLAS_DYNAMIC_ARCH)
            set(OPENBLAS_DYNAMIC_ARCH ON)
        endif()
        if(NOT DEFINED OPENBLAS_TARGET)
            set(OPENBLAS_TARGET HASWELL)
        endif()
        install_package(OpenBLAS VERSION 0.3.17
                DEPENDS gfortran::gfortran Threads::Threads
                CMAKE_ARGS
                -DDYNAMIC_ARCH:BOOL=${OPENBLAS_DYNAMIC_ARCH}
                -DTARGET:BOOL=${OPENBLAS_TARGET}
                -DUSE_THREAD:BOOL=ON
                -DBUILD_RELAPACK:BOOL=OFF
                )
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        # Fix for OpenBLAS 0.3.9, which otherwise includes <complex> inside of an extern "C" scope.
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
        #For convenience, define these targes
        if (NOT TARGET BLAS::BLAS)
            add_library(BLAS::BLAS INTERFACE IMPORTED)
            target_link_libraries(BLAS::BLAS INTERFACE OpenBLAS::OpenBLAS)
        endif ()
        if (NOT TARGET LAPACK::LAPACK)
            add_library(LAPACK::LAPACK INTERFACE IMPORTED)
            target_link_libraries(LAPACK::LAPACK INTERFACE OpenBLAS::OpenBLAS)
        endif ()
    endif ()


    # Starting from here there should definitely be blas library that includes lapacke
    # Lapacke is needed by arpack++, included in MKL or OpenBLAS
    find_package(Lapacke REQUIRED)

    # cxxopts for parsing cli arguments
    install_package(cxxopts VERSION 2.2.0)

    # Eigen3 numerical library (needed by ceres and h5pp)
    install_package(Eigen3 VERSION 3.4 TARGET_NAME Eigen3::Eigen)
    # h5pp for writing to file binary in format
    install_package(h5pp VERSION 1.9.0 CMAKE_ARGS ${h5pp_ARGS})
    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    install_package(arpack-ng VERSION 3.8.0
            TARGET_NAME ARPACK::ARPACK
            DEPENDS BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)

    # C++ frontend for arpack-ng. Custom find module.
    install_package(arpack++
            TARGET_NAME arpack::arpack++
            DEPENDS ARPACK::ARPACK lapacke::lapacke gfortran::gfortran
            MODULE CHECK)

    # Google Flags library needed by ceres-solver
    install_package(gflags VERSION 2.2.2 COMPONENTS ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS})

    # Google logging library needed by ceres-solver
    install_package(glog VERSION 0.5 CMAKE_ARGS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR} CHECK)

    # ceres-solver (for L-BFGS routine)
    install_package(Ceres VERSION 2.0
            TARGET_NAME Ceres::ceres
            DEPENDS gflags glog::glog
            CMAKE_ARGS ${Ceres_ARGS}
            CHECK
            QUIET)


    # Configure Eigen

    if(TARGET Eigen3::Eigen AND DMRG_ENABLE_THREADS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    endif()

    if(TARGET Eigen3::Eigen)
        get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})


        if(TARGET mkl::mkl)
            message(STATUS "Eigen3 will use MKL")
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
            target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
        elseif(TARGET OpenBLAS::OpenBLAS)
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

endif ()
