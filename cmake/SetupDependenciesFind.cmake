if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    # Add search directories and flags for the CMake find_* tools
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    if(DMRG_PACKAGE_MANAGER MATCHES "cmake")
        # We set variables here that allows us to find packages exclusively with CMAKE_PREFIX_PATH
        list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX})
        list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
        set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE STRING "" FORCE)
        set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
        # Flags that can be used directly on find_package
        # Enumerated according to the cmake manual for find_package
        # These lets us ignore system packages when pkg manager matches "cmake"
        set(N5 NO_SYSTEM_ENVIRONMENT_PATH) #5
        set(N6 NO_CMAKE_PACKAGE_REGISTRY) #6
        set(N7 NO_CMAKE_SYSTEM_PATH) #7
        set(N8 NO_CMAKE_SYSTEM_PACKAGE_REGISTRY) #8
    endif()
    if(NOT BUILD_SHARED_LIBS)
        set(GFLAGS_COMPONENTS COMPONENTS)
        set(GFLAS_ITEMS nothreads_static)
    endif()

    if(DMRG_ENABLE_OPENMP)
        find_package(OpenMP COMPONENTS CXX REQUIRED)
        set(mkl_thread gnu_thread)
    else()
        set(mkl_thread sequential)
    endif()

    find_package(Fortran REQUIRED)

    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS gf ${mkl_thread} lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    else()
        find_package(OpenBLAS 0.3.8 ${N5} ${N6} ${N7} ${N8} ${REQUIRED}) # If MKL is not on openblas will be used instead. Includes lapack.
    endif()
    find_package(Lapacke REQUIRED)                                      # Lapacke needed by arpack++, included in MKL or OpenBLAS

    find_package(Eigen3 3.3.7 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})      # Eigen3 numerical library (needed by ceres and h5pp)
    find_package(h5pp 1.9.0 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})        # h5pp for writing to file binary in format
    find_package(arpack-ng 3.8.0 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})   # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.

    # Arpack needs to link to extra libraries before moving on
    if(arpack-ng_FOUND)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()

    find_package(arpack++ ${REQUIRED})                                  # C++ frontend for arpack-ng. Custom find module.
    find_package(gflags 2.2.2 ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS}
                              ${N5} ${N6} ${N7} ${N8} ${REQUIRED})      # Google Flags library needed by ceres-solver
    find_package(glog 0.4 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})          # Google logging library needed by ceres-solver
    find_package(Ceres 2.0 PATH_SUFFIXES ceres ceres/lib
                              ${N5} ${N6} ${N7} ${N8} ${REQUIRED})      # ceres-solver (for L-BFGS routine)





endif()