cmake_minimum_required(VERSION 3.15)
project(external-glog)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)
include(ExternalProject)
set(cxxflags "-w")
ExternalProject_Add(external_glog
        URL https://github.com/google/glog/archive/refs/tags/v0.6.0.tar.gz
        URL_MD5 c98a6068bc9b8ad9cebaca625ca73aa2
        BUILD_ALWAYS TRUE
        PREFIX ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
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
        # glog flags
        -Dgflags_ROOT:PATH=${gflags_ROOT}
        -DWITH_GTEST:BOOL=FALSE
        -DWITH_PKGCONFIG:BOOL=FALSE
        -DWITH_SYMBOLIZE:BOOL=FALSE
        -DWITH_THREADS:BOOL=TRUE
        -DWITH_TLS:BOOL=TRUE
        -DWITH_UNWIND:BOOL=FALSE
        )