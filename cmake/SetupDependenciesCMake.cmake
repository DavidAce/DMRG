
if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    include(cmake/InstallPackage.cmake)

    # Only give one chance to find the package if pkg manager is "find"
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()

    # We set variables here that allows us to find packages with CMAKE_PREFIX_PATH
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

    # Set CMake build options
    list(APPEND OpenBLAS_CMAKE_OPTIONS -DTARGET:STRING=${OPENBLAS_MARCH})
    list(APPEND OpenBLAS_CMAKE_OPTIONS -DUSE_THREAD:BOOL=1)
    list(APPEND OpenBLAS_CMAKE_OPTIONS -DBUILD_RELAPACK:BOOL=OFF)

    list(APPEND h5pp_CMAKE_OPTIONS -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND h5pp_CMAKE_OPTIONS -DH5PP_PACKAGE_MANAGER:STRING=cmake)
    list(APPEND h5pp_CMAKE_OPTIONS -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})

    list(APPEND glog_CMAKE_OPTIONS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})

    list(APPEND Ceres_CMAKE_OPTIONS -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND Ceres_CMAKE_OPTIONS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
    list(APPEND Ceres_CMAKE_OPTIONS -Dglog_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})

    if(NOT BUILD_SHARED_LIBS)
        set(GFLAGS_COMPONENTS COMPONENTS)
        set(GFLAS_ITEMS nothreads_static)
    endif()




    # Find packages or install if missing

    if(DMRG_ENABLE_OPENMP)
        find_package(OpenMP COMPONENTS CXX REQUIRED)
        set(mkl_thread gnu_thread)
    else()
        set(mkl_thread sequential)
    endif()

    find_package(Fortran REQUIRED)

    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS blas lapack gf ${mkl_thread} lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
    else()
        find_package(OpenBLAS 0.3.8 ${N5} ${N6} ${N7} ${N8} ${REQUIRED}) # If MKL is not on openblas will be used instead. Includes lapack.
    endif()

    if(NOT MKL_FOUND)
        install_package(OpenBLAS "${DMRG_DEPS_INSTALL_DIR}" "${OpenBLAS_CMAKE_OPTIONS}")
        find_package(OpenBLAS 0.3.8 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()
    if(TARGET OpenBLAS::OpenBLAS)
        target_link_libraries(OpenBLAS::OpenBLAS INTERFACE gfortran::gfortran Threads::Threads)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        # Fix for OpenBLAS 0.3.9, which otherwise includes <complex> inside of an extern "C" scope.
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
        #For convenience, define these targes
        if(NOT TARGET BLAS::BLAS)
            add_library(BLAS::BLAS                  INTERFACE IMPORTED)
            target_link_libraries(BLAS::BLAS        INTERFACE OpenBLAS::OpenBLAS)
        endif()
        if(NOT TARGET LAPACK::LAPACK)
            add_library(LAPACK::LAPACK              INTERFACE IMPORTED)
            target_link_libraries(LAPACK::LAPACK    INTERFACE OpenBLAS::OpenBLAS)
        endif()
    endif()


    # Starting from here there should definitely be blas library that includes lapacke
    # Lapacke is needed by arpack++, included in MKL or OpenBLAS
    find_package(Lapacke REQUIRED)

    # Eigen3 numerical library (needed by ceres and h5pp)
    find_package(Eigen3 3.3.7 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    install_package(Eigen3 "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(Eigen3 3.3.7 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    # h5pp for writing to file binary in format
    find_package(h5pp 1.9.0 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    install_package(h5pp "${DMRG_DEPS_INSTALL_DIR}" "${h5pp_CMAKE_OPTIONS}")
    find_package(h5pp 1.9.1 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
    find_package(arpack-ng 3.8.0 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    if(NOT arpack-ng_FOUND)
        install_package(arpack-ng "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(arpack-ng 3.8.0 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()
    target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)

    # C++ frontend for arpack-ng. Custom find module.
    find_package(arpack++ ${REQUIRED})
    install_package(arpack++ "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(arpack++ REQUIRED)

    # Google Flags library needed by ceres-solver
    find_package(gflags 2.2.2 ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS} ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    install_package(gflags "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(gflags 2.2.2 ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS} HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    # Google logging library needed by ceres-solver
    find_package(glog 0.4 ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    install_package(glog "${DMRG_DEPS_INSTALL_DIR}" "${glog_CMAKE_OPTIONS}")
    find_package(glog 0.4 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    include(cmake/CheckGlogCompiles.cmake)
    check_glog_compiles("glog::glog" "" "" "" "")

    # ceres-solver (for L-BFGS routine)
    find_package(Ceres 2.0 PATH_SUFFIXES ceres ceres/lib ${N5} ${N6} ${N7} ${N8} ${REQUIRED})
    if(NOT Ceres_FOUND)
        install_package(Ceres "${DMRG_DEPS_INSTALL_DIR}" "${Ceres_CMAKE_OPTIONS}" )
        find_package(Ceres 2.0 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
        include(cmake/CheckCeresCompiles.cmake)
        check_ceres_compiles("Ceres::ceres" "" "" "" "")
    endif()

endif()
