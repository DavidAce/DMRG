cmake_minimum_required(VERSION 3.15)
project(external-Ceres)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)
option(CMAKE_VERBOSE_MAKEFILE "Verbosity" OFF)
option(GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION        ""  TRUE)
option(EIGEN_PREFER_EXPORTED_EIGEN_CMAKE_CONFIGURATION      ""  TRUE)
option(GFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION    ""  TRUE)


include(ExternalProject)
ExternalProject_Add(external_Ceres
        URL  https://github.com/ceres-solver/ceres-solver/archive/2.0.0.tar.gz
        URL_MD5 94246057ac520313e3b582c45a30db6e
        BUILD_ALWAYS TRUE # <- True will detect changes in dependencies of ceres
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_SKIP_BUILD_RPATH:BOOL=OFF
        -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
        -DCMAKE_INSTALL_MESSAGE=LAZY
        # Ceres flags
        -DEigen3_ROOT:PATH=${Eigen3_ROOT}
        -Dgflags_ROOT:PATH=${gflags_ROOT}
        -Dglog_ROOT:PATH=${glog_ROOT}
        -DCERES_THREADING_MODEL=NO_THREADS
        -DBUILD_TESTING:BOOL=OFF
        -DBUILD_EXAMPLES:BOOL=OFF
        -DLAPACK:BOOL=OFF
        -DEIGENSPARSE:BOOL=ON
        -DGFLAGS:BOOL=ON
        -DSUITESPARSE:BOOL=OFF
        -DCXSPARSE:BOOL=OFF
        -DSCHUR_SPECIALIZATIONS:BOOL=OFF
        -DCUSTOM_BLAS:BOOL=ON
        )