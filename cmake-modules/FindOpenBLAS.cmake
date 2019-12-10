function(find_OpenBLAS)
    # It seems that libopenblas from "apt" includes lapack in newer versions of Ubuntu. Trusty LTS 14.04 does not
    # bundle OpenBLAS with lapack. Therefore we test if both lapack and blas are present to distinghuish these cases.
    # Otherwise, arpack-ng will complain later about undefined references.
    set(BLA_VENDOR OpenBLAS)
    set(BLAS_VERBOSE OFF)
    if(NOT ${BUILD_SHARED_LIBS})
        set(BLA_STATIC ON)
    else()
        set(BLA_STATIC OFF)
    endif()


    if(NOT BLAS_FOUND)
        message(STATUS "Searching for OpenBLAS config")
        find_package(OpenBLAS 0.3
                HINTS
                $ENV{CONDA_PREFIX}
                ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                $ENV{EBROOTOPENBLAS}
                PATHS
                ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                $ENV{EBROOTOPENBLAS}
                $ENV{CONDA_PREFIX}
                $ENV{BLAS_DIR}/lib
                /usr/lib/x86_64-linux-gnu/
                /usr/lib/
                PATH_SUFFIXES
                lib OpenBLAS/lib
                )

        if(OpenBLAS_FOUND AND OpenBLAS_LIBRARIES AND NOT OpenBLAS_LIBRARY)
            get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
            get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
            if(EXISTS  ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
                set(OpenBLAS_LIBRARY ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
            else()
                set(STATUS "Found OpenBLAS library with wrong suffix (i.e. not static/shared as requested")
            endif()
        endif()

        if(OpenBLAS_LIBRARY AND OpenBLAS_FOUND AND OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIRS)
            set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
            set(OpenBLAS_FOUND TRUE)
            set(BLAS_FOUND TRUE)
        endif()
    endif()




    if(NOT OpenBLAS_FOUND)
        message(STATUS "Searching for OpenBLAS in library file")
        find_library(OpenBLAS_LIBRARY
                NAMES libopenblas${CUSTOM_SUFFIX}
                HINTS
                $ENV{CONDA_PREFIX}
                ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                $ENV{EBROOTOPENBLAS}
                PATHS
                ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                $ENV{EBROOTOPENBLAS}
                $ENV{CONDA_PREFIX}
                $ENV{BLAS_DIR}/lib
                /usr/lib/x86_64-linux-gnu/
                /usr/lib/
                PATH_SUFFIXES
                lib OpenBLAS/lib

                )
        find_path(OpenBLAS_INCLUDE_DIRS
                NAMES openblas_config.h
                PATHS
                ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                $ENV{EBROOTOPENBLAS}/include
                $ENV{CONDA_PREFIX}
                $ENV{BLAS_DIR}/include
                /usr/include
                /usr/include/x86_64-linux-gnu
                PATH_SUFFIXES
                include OpenBLAS/include blas/include
                )
        if (OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIRS)
            # External libs are probably multithreaded
            set(OpenBLAS_MULTITHREADED 1 )
            set(BLAS_FOUND TRUE)
            set(OpenBLAS_FOUND TRUE)
        endif()
    endif()
    if(OpenBLAS_FOUND)
        message(STATUS "OpenBLAS FOUND IN SYSTEM: ${OpenBLAS_LIBRARY}")
        message(STATUS "                          ${OpenBLAS_INCLUDE_DIRS}")
        #For convenience, define these variables
        add_library(OpenBLAS   INTERFACE IMPORTED)
        target_link_libraries(OpenBLAS INTERFACE ${OpenBLAS_LIBRARY} gfortran Threads::Threads)
        target_include_directories(OpenBLAS INTERFACE ${OpenBLAS_INCLUDE_DIRS})
        target_compile_definitions(OpenBLAS INTERFACE OpenBLAS_AVAILABLE)
        #For convenience, define these variables
        add_library(blas   INTERFACE IMPORTED)
        add_library(lapack INTERFACE IMPORTED)
        target_link_libraries(blas   INTERFACE OpenBLAS)
        target_link_libraries(lapack INTERFACE OpenBLAS)
    endif()
endfunction()