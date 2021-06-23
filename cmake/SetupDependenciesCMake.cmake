
# Make sure to include cmake/SetupDependenciesFind.cmake before this


if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    include(cmake/InstallPackage.cmake)

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


    # Install missing packages and try finding them again
    if(NOT MKL_FOUND)
        install_package(OpenBLAS "${DMRG_DEPS_INSTALL_DIR}" "${OpenBLAS_CMAKE_OPTIONS}")
        find_package(OpenBLAS 0.3.8 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()

    install_package(Eigen3 "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(Eigen3 3.3.7 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    install_package(h5pp "${DMRG_DEPS_INSTALL_DIR}" "${h5pp_CMAKE_OPTIONS}")
    find_package(h5pp 1.9.1 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    if(NOT arpack-ng_FOUND)
        install_package(arpack-ng "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(arpack-ng 3.8.0 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()



    install_package(arpack++ "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(arpack++ REQUIRED)

    install_package(gflags "${DMRG_DEPS_INSTALL_DIR}" "")
    find_package(gflags ${GFLAGS_COMPONENTS} ${GFLAGS_ITEMS} HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    install_package(glog "${DMRG_DEPS_INSTALL_DIR}" "${glog_CMAKE_OPTIONS}")
    find_package(glog 0.4 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)

    install_package(Ceres "${DMRG_DEPS_INSTALL_DIR}" "${Ceres_CMAKE_OPTIONS}" )
    find_package(Ceres HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)


endif()
