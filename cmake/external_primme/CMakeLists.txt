cmake_minimum_required(VERSION 3.15)
project(external-primme)
include(../GetNumThreads.cmake)
get_num_threads(NUM_THREADS)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(ExternalProject)

set(cflags "-O3 -DNDEBUG -std=c17 -fPIC -g -fno-strict-aliasing -fno-omit-frame-pointer ${CMAKE_CXX_FLAGS_RELEASE} -DPRIMME_BLASINT_SIZE=32 -DPRIMME_INT_SIZE=0")
set(cxxflags "-O3 -DNDEBUG -std=c++17 -fPIC -g -fno-strict-aliasing -fno-omit-frame-pointer ${CMAKE_C_FLAGS_RELEASE} -DPRIMME_BLASINT_SIZE=32 -DPRIMME_INT_SIZE=0")

if(BUILD_SHARED_LIBS)
    set(type solib)
else()
    set(type lib)
endif()

externalproject_add(external_primme
                    URL https://github.com/primme/primme/archive/refs/tags/v3.2.tar.gz
                    URL_MD5 b26968d0ea8aa2e6feefc89f3c863062
                    PREFIX ${CMAKE_BINARY_DIR}
                    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}

                    #Give Ninja access to the terminal.
                    USES_TERMINAL_DOWNLOAD TRUE
                    USES_TERMINAL_UPDATE TRUE
                    USES_TERMINAL_CONFIGURE TRUE
                    USES_TERMINAL_BUILD TRUE
                    USES_TERMINAL_INSTALL TRUE
                    USES_TERMINAL_TEST TRUE

                    BUILD_ALWAYS TRUE
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND
                    COMMAND cd <SOURCE_DIR>
                    COMMAND make VERBOSE=1 ${type} -j${NUM_THREADS} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} CXXFLAGS=${cxxflags} CFLAGS=${cflags}

                    INSTALL_COMMAND ""
                    )

if(BUILD_SHARED_LIBS)
    externalproject_add_step(
            external_primme
            install_shared_lib
            DEPENDEES build
            COMMENT "Installing shared lib"
            COMMAND cd <SOURCE_DIR>
            COMMAND make install includedir=<INSTALL_DIR>/include/primme libdir=<INSTALL_DIR>/lib
    )
else()
    externalproject_add_step(
            external_primme
            install_static_lib
            DEPENDEES build
            COMMENT "Installing static lib"
            COMMAND cmake -E copy <SOURCE_DIR>/lib/libprimme${CMAKE_STATIC_LIBRARY_SUFFIX} <INSTALL_DIR>/lib
    )
endif()
