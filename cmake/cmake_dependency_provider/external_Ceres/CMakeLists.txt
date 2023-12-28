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

set(cxxflags "-w")

include(ExternalProject)
ExternalProject_Add(external_Ceres
        #URL https://github.com/ceres-solver/ceres-solver/archive/refs/tags/2.1.0.tar.gz
        URL https://github.com/ceres-solver/ceres-solver/archive/refs/tags/2.2.0.tar.gz
        #URL_MD5 4ab48f439d63a09a075d4fcdcddc7667 #2.1.0
        URL_MD5 b360c795a75c3e99d67a2a60f2004466  #2.2.0
        BUILD_ALWAYS TRUE # <- True will detect changes in dependencies of ceres
        PREFIX ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE

        #Give Ninja access to the terminal.
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_UPDATE TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
        USES_TERMINAL_TEST TRUE
        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
        #-DCMAKE_SKIP_BUILD_RPATH:BOOL=OFF
        #-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
        -DCMAKE_INSTALL_MESSAGE=LAZY
        -DCMAKE_CXX_FLAGS=${cxxflags}
        -DBUILD_SHARED_LIBS:BOOL=OFF
        # Ceres flags
        -DMINIGLOG:BOOL=TRUE
        -DBUILD_TESTING:BOOL=OFF
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_DOCUMENTATION:BOOL=OFF
        -DBUILD_BENCHMARKS:BOOL=OFF
        -DLAPACK:BOOL=OFF
        -DEIGENSPARSE:BOOL=ON
        -DEIGENMETIS:BOOL=OFF
        -DGFLAGS:BOOL=OFF
        -DSUITESPARSE:BOOL=OFF
        -DCXSPARSE:BOOL=OFF
        -DSCHUR_SPECIALIZATIONS:BOOL=OFF
        -DCUSTOM_BLAS:BOOL=ON
        -DUSE_CUDA:BOOL=OFF
        -DEigen3_ROOT:PATH=<INSTALL_DIR>

        # These are only needed if we enable glog and gflags
        -Dgflags_ROOT:PATH=<INSTALL_DIR>
        -Dglog_ROOT:PATH=<INSTALL_DIR>
        -DGLOG_FOUND:BOOL=TRUE # To fix a bug with Ceres looking checking GLOG_FOUND instead of Glog_FOUND
        )