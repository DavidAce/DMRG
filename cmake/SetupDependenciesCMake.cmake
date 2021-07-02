
if (DMRG_PACKAGE_MANAGER STREQUAL "cmake")
    include(cmake/InstallPackage.cmake)

    # Set CMake build options
    if(OPENBLAS_TARGET)
        list(APPEND OpenBLAS_ARGS -DTARGET:STRING=${OPENBLAS_MARCH})
    endif()
    list(APPEND OpenBLAS_ARGS -DUSE_THREAD:BOOL=ON)
    list(APPEND OpenBLAS_ARGS -DBUILD_RELAPACK:BOOL=OFF)

    list(APPEND h5pp_ARGS -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND h5pp_ARGS -DH5PP_PACKAGE_MANAGER:STRING=cmake)
    list(APPEND h5pp_ARGS -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})

    list(APPEND glog_ARGS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})

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
        install_package(OpenBLAS VERSION 0.3.8
                CMAKE_ARGS ${OpenBLAS_ARGS}
                DEPENDS gfortran::gfortran Threads::Threads)
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

    # Eigen3 numerical library (needed by ceres and h5pp)
    install_package(Eigen3 VERSION 3.3)
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
    install_package(glog VERSION 0.4 CMAKE_ARGS ${glog_ARGS} CHECK)

    # ceres-solver (for L-BFGS routine)
    install_package(Ceres VERSION 2.0
            TARGET_NAME Ceres::ceres
            DEPENDS gflags glog::glog
            CMAKE_ARGS ${Ceres_ARGS}
            CHECK)
endif ()
