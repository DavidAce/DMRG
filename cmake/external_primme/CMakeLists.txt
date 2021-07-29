cmake_minimum_required(VERSION 3.15)
project(external-primme)
include(../GetNumThreads.cmake)
get_num_threads(NUM_THREADS)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(ExternalProject)
#set(cc "$ENV{CXX}")
#set(cflags "-O3 -DNDEBUG -fPIC")
set(cflags "-O3 -std=c17 -fPIC -g -fno-strict-aliasing -fno-omit-frame-pointer -DPRIMME_BLASINT_SIZE=32 -DPRIMME_INT_SIZE=0")
set(cxxflags "-O3 -std=c++17 -fPIC -g -fno-strict-aliasing -fno-omit-frame-pointer -DPRIMME_BLASINT_SIZE=32 -DPRIMME_INT_SIZE=0")
if(BUILD_SHARED_LIBS)
    set(type solib)
else()
    set(type lib)
endif()
ExternalProject_Add(external_primme
        URL      https://github.com/primme/primme/archive/refs/tags/v3.2.tar.gz
        URL_MD5  b26968d0ea8aa2e6feefc89f3c863062
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS TRUE
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
            cd <SOURCE_DIR> &&
            pwd &&
            make ${type} -j${NUM_THREADS}  CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} CXXFLAGS=${cxxflags} CFLAGS=${cflags}

        INSTALL_COMMAND
            cmake -E copy_directory <SOURCE_DIR>/lib <INSTALL_DIR>/lib &&
            cmake -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include/primme
        )