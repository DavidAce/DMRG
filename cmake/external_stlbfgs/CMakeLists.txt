cmake_minimum_required(VERSION 3.15)
project(external-stlbfgs)
if (INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif ()
include(GNUInstallDirs)
include(ExternalProject)
set(cxxflags "-w")

if (BUILD_SHARED_LIBS)
    set(STLBFGS_LIBSUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else ()
    set(STLBFGS_LIBSUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif ()

ExternalProject_Add(external_stlbfgs
        GIT_REPOSITORY https://github.com/ssloy/stlbfgs.git
        GIT_TAG 0742ff7c4c0d76656a2a9c6723b7b13e5a834aeb
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

        INSTALL_COMMAND
        cmake -E make_directory <INSTALL_DIR>/lib &&
        cmake -E make_directory <INSTALL_DIR>/include/stlbfgs &&
        cmake -E copy_if_different <BINARY_DIR>/libstlbfgs${STLBFGS_LIBSUFFIX} <INSTALL_DIR>/lib &&
        cmake -E copy_if_different <SOURCE_DIR>/linesearch.h <INSTALL_DIR>/include/stlbfgs/linesearch.h &&
        cmake -E copy_if_different <SOURCE_DIR>/stlbfgs.h <INSTALL_DIR>/include/stlbfgs/stlbfgs.h
        )