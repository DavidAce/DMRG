
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
            set(OPENBLAS_DYNAMIC_ARCH OFF)
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
            message(STATUS "Defined target BLAS::BLAS -> OpenBLAS::OpenBLAS")
        endif ()
        if (NOT TARGET LAPACK::LAPACK)
            add_library(LAPACK::LAPACK INTERFACE IMPORTED)
            target_link_libraries(LAPACK::LAPACK INTERFACE OpenBLAS::OpenBLAS)
            message(STATUS "Defined target LAPACK::LAPACK -> OpenBLAS::OpenBLAS")
        endif ()
    endif ()

    # Starting from here there should definitely be blas library that includes lapacke
    # Lapacke is needed by arpack++, included in MKL or OpenBLAS
    find_package(Lapacke REQUIRED)

    install_package(pcg-cpp)

    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    install_package(arpack-ng VERSION 3.8.0
            TARGET_NAME ARPACK::ARPACK
            DEPENDS BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)

    # C++ frontend for arpack-ng. Custom find module.
    install_package(arpack++
            TARGET_NAME arpack++::arpack++
            DEPENDS ARPACK::ARPACK lapacke::lapacke gfortran::gfortran
            MODULE CHECK)

    # Eigen3 numerical library (needed by ceres and h5pp)
    install_package(Eigen3 VERSION 3.4 TARGET_NAME Eigen3::Eigen)

    # cli11 for parsing cli arguments
    install_package(cli11 VERSION 2.1.1 TARGET_NAME CLI11::CLI11 FIND_NAME CLI11)


    # h5pp for writing to file binary in format
    install_package(h5pp VERSION 1.10.0 CMAKE_ARGS ${h5pp_ARGS})

    # Backward for printing pretty stack traces
    install_package(Backward MODULE)

    # Google Flags library needed by ceres-solver
    install_package(gflags VERSION 2.2.2 COMPONENTS ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS})

    # Google logging library needed by ceres-solver
    install_package(glog VERSION 0.6
            #            DEPENDS libunwind::libunwind
            CMAKE_ARGS
            -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR}
            #                -DUnwind_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR}
            CHECK)

    # ceres-solver (for L-BFGS routine)
    install_package(Ceres VERSION 2.1.0
            TARGET_NAME Ceres::ceres
            DEPENDS gflags glog::glog
            CMAKE_ARGS ${Ceres_ARGS}
            CHECK
            QUIET)

    target_link_libraries(dmrg-deps INTERFACE
            CLI11::CLI11
            pcg-cpp::pcg-cpp
            h5pp::h5pp
            arpack++::arpack++
            Ceres::ceres
            BLAS::BLAS
            Backward::Backward
            )
endif ()
