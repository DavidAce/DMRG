cmake_minimum_required(VERSION 3.15)
project(external-gflags)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)
include(ExternalProject)
ExternalProject_Add(external_gflags
        URL https://github.com/gflags/gflags/archive/v2.2.2.tar.gz
        URL_MD5 1a865b93bacfa963201af3f75b7bd64c
        BUILD_ALWAYS TRUE
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        UPDATE_COMMAND ""
        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_SKIP_BUILD_RPATH:BOOL=OFF
        -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
        -DCMAKE_INSTALL_MESSAGE=LAZY
        # gflags flags
        -DGFLAGS_BUILD_SHARED_LIBS:BOOL=ON
        -DGFLAGS_BUILD_STATIC_LIBS:BOOL=ON
        )