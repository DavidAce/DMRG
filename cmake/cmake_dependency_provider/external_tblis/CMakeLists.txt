cmake_minimum_required(VERSION 3.15)
project(external-tblis)
cmake_host_system_information(RESULT NUM_THREADS QUERY NUMBER_OF_PHYSICAL_CORES)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(ExternalProject)

set(cc ${CMAKE_C_COMPILER})
set(cxx ${CMAKE_CXX_COMPILER})
set(cflags "-std=c17 -fPIC ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
set(cxxflags "-std=c++17 -fPIC ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
set(ldflags "${CMAKE_EXE_LINKER_FLAGS}")
set(libs "-lm -ldl")
if(TBLIS_THREAD_MODEL)
    set(thread_model "--enable-thread-model=${TBLIS_THREAD_MODEL}")
endif()
externalproject_add(external_tblis
                    GIT_REPOSITORY https://github.com/devinamatthews/tblis.git
                    GIT_TAG 7d979b6624071879383377611cb2fa483b6a3104 # master on Apr 22
                    #         GIT_TAG  3845b1744c186ca3be45f9c9a69bcf527ac754b6 # develop on May 15, 2021
                    GIT_PROGRESS TRUE
                    #        PATCH_COMMAND patch -p1 -i ${PROJECT_SOURCE_DIR}/tblis-develop-3845b17.patch
                    PATCH_COMMAND patch -p1 -i ${PROJECT_SOURCE_DIR}/tblis-config.patch || true # For 7d979b66

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
                    BUILD_ALWAYS TRUE
                    BUILD_IN_SOURCE TRUE
                    CONFIGURE_COMMAND
                    COMMAND export CC=${cc}
                    COMMAND export CXX=${cxx}
                    COMMAND export LIBS=${libs}
                    COMMAND export CFLAGS=${cflags}
                    COMMAND export CXXFLAGS=${cxxflags}
                    COMMAND export LDFLAGS=${ldflags}
                    COMMAND ./configure --enable-config=haswell,skx,zen ${thread_model} --prefix=<INSTALL_DIR>

                    BUILD_COMMAND
                    COMMAND make -j${NUM_THREADS}

                    INSTALL_COMMAND
                    make install
                    )