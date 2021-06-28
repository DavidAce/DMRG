#cmake_minimum_required(VERSION 3.15)
#project(external-glog)
#if(INIT_CACHE_FILE)
#    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
#endif()
#include(GNUInstallDirs)
#include(ExternalProject)
#ExternalProject_Add(external_glog
#        URL https://github.com/google/glog/archive/v0.4.0.tar.gz
#        URL_MD5 0daea8785e6df922d7887755c3d100d0
#        BUILD_ALWAYS TRUE
#        PREFIX          ${CMAKE_BINARY_DIR}
#        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
#        CMAKE_ARGS
#        ${INIT_CACHE_ARG}
#        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
#        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
#        -DCMAKE_SKIP_BUILD_RPATH:BOOL=OFF
#        -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
#        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
#        -DCMAKE_INSTALL_MESSAGE=LAZY
#        # glog flags
#        -Dgflags_ROOT:PATH=${gflags_ROOT}
#        -DBUILD_TESTING:BOOL=TRUE
#        -DWITH_THREADS:BOOL=TRUE
#        -DWITH_TLS:BOOL=TRUE
#        )



cmake_minimum_required(VERSION 3.15)
project(external-glog)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)
include(ExternalProject)
ExternalProject_Add(external_glog
        URL https://github.com/google/glog/archive/v0.5.0.tar.gz
        URL_MD5 2368e3e0a95cce8b5b35a133271b480f
        BUILD_ALWAYS TRUE
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
        # glog flags
        -Dgflags_ROOT:PATH=${gflags_ROOT}
        -DWITH_GTEST:BOOL=FALSE
        -DWITH_PKGCONFIG:BOOL=FALSE
        -DWITH_SYMBOLIZE:BOOL=FALSE
        -DWITH_THREADS:BOOL=TRUE
        -DWITH_TLS:BOOL=TRUE
        -DWITH_UNWIND:BOOL=FALSE
        )