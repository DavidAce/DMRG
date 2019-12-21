function(find_OpenBLAS)
    if(NOT TARGET openblas::openblas)
        message(STATUS "Searching for OpenBLAS config")
        find_package(OpenBLAS 0.3
                HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS}
                PATHS
                    $ENV{EBROOTBLAS}
                    $ENV{BLAS_DIR}
                    $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                    lib OpenBLAS/lib OpenBLAS openblas
                )

        if(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            get_filename_component(OpenBLAS_LIBRARIES_EXT ${OpenBLAS_LIBRARIES} EXT)
            if(NOT "${OpenBLAS_LIBRARIES_EXT}" MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}" AND NOT BUILD_SHARED_LIBS)
                get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
                get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
                if(EXISTS  ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX} )
                    message(STATUS "Replaced OpenBLAS extension from ${OpenBLAS_LIBRARIES_EXT} to ${CMAKE_STATIC_LIBRARY_SUFFIX}")
                    set(OpenBLAS_LIBRARIES ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX})
                else()
                    message(STATUS "Could not find OpenBLAS library with correct suffix (i.e. static/shared as requested): ${OpenBLAS_LIBRARIES}")
                endif()
            endif()
            if(NOT TARGET openblas::openblas)
                add_library(openblas::openblas ${LINK_TYPE} IMPORTED)
                set_target_properties(openblas::openblas PROPERTIES
                        IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                        INTERFACE_COMPILE_DEFINITIONS "OpenBLAS_AVAILABLE")
                target_include_directories(openblas::openblas SYSTEM INTERFACE ${OpenBLAS_INCLUDE_DIRS})

            endif()
        endif()
    endif()




    if(NOT TARGET openblas::openblas)
        message(STATUS "Searching for OpenBLAS in system")
        find_library(OpenBLAS_LIBRARIES
                NAMES openblas
                HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS}
                PATHS
                $ENV{EBROOTBLAS}
                ${CONDA_HINTS}
                $ENV{BLAS_DIR}
                $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                    lib openblas/lib OpenBLAS/lib openblas OpenBLAS

                )
        find_path(OpenBLAS_INCLUDE_DIRS
                NAMES openblas_config.h
                HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS}
                PATHS
                    $ENV{EBROOTBLAS}
                    $ENV{BLAS_DIR}
                    $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                    include openblas openblas/include OpenBLAS OpenBLAS/include blas/include
                )
        if (OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            add_library(openblas::openblas STATIC IMPORTED)
            set_target_properties(openblas::openblas PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                    INTERFACE_COMPILE_DEFINITIONS "OpenBLAS_AVAILABLE")
            target_include_directories(openblas::openblas SYSTEM INTERFACE ${OpenBLAS_INCLUDE_DIRS})


        endif()
    endif()

    if(TARGET openblas::openblas)
        message(STATUS "OpenBLAS found: ${OpenBLAS_LIBRARIES}")
        message(STATUS "                ${OpenBLAS_INCLUDE_DIRS}")

        target_link_libraries(openblas::openblas INTERFACE gfortran::gfortran pthread)
        #For convenience, define these targes
        add_library(blas::blas       INTERFACE IMPORTED)
        add_library(lapack::lapack   INTERFACE IMPORTED)
        target_link_libraries(blas::blas        INTERFACE openblas::openblas)
        target_link_libraries(lapack::lapack    INTERFACE openblas::openblas)
    else()
        message(STATUS "Could not find OpenBLAS")
    endif()
endfunction()