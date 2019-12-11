function(find_OpenBLAS)
    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS config")
        find_package(OpenBLAS 0.3
                HINTS
                    ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                    $ENV{CONDA_PREFIX}
                    $ENV{EBROOTOPENBLAS}
                PATHS
                    ${CMAKE_INSTALL_PREFIX}/OpenBLAS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{CONDA_PREFIX}
                    $ENV{BLAS_DIR}/lib
                PATH_SUFFIXES
                    lib OpenBLAS/lib
                )

        if(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
            get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
            if(EXISTS  ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
                set(OpenBLAS_LIBRARIES ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
                set(OpenBLAS_FOUND TRUE)
                set(BLAS_FOUND TRUE)
                add_library(OpenBLAS INTERFACE IMPORTED)
                target_link_libraries(OpenBLAS INTERFACE ${OpenBLAS_LIBRARIES})
                target_include_directories(OpenBLAS INTERFACE ${OpenBLAS_INCLUDE_DIRS})
            else()
                message(STATUS "Found OpenBLAS library with wrong suffix (i.e. not static/shared as requested")
            endif()
        endif()
    endif()




    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS lib in system")
        find_library(OpenBLAS_LIBRARIES
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
        if (OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            add_library(OpenBLAS INTERFACE IMPORTED)
            target_link_libraries(OpenBLAS INTERFACE ${OpenBLAS_LIBRARIES})
            target_include_directories(OpenBLAS INTERFACE ${OpenBLAS_INCLUDE_DIRS})
        endif()
    endif()

    if(TARGET OpenBLAS)
        message(STATUS "OpenBLAS found: ${OpenBLAS_LIBRARIES}")
        message(STATUS "                ${OpenBLAS_INCLUDE_DIRS}")
        target_link_libraries(OpenBLAS INTERFACE gfortran Threads::Threads)
        target_compile_definitions(OpenBLAS INTERFACE OpenBLAS_AVAILABLE)
        #For convenience, define these variables
        add_library(blas   INTERFACE IMPORTED)
        add_library(lapack INTERFACE IMPORTED)
        target_link_libraries(blas   INTERFACE OpenBLAS)
        target_link_libraries(lapack INTERFACE OpenBLAS)
    else()
        message(STATUS "Could not find OpenBLAS")
    endif()
endfunction()