function(find_OpenBLAS)
    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS config")
        message(STATUS "DIRECTORY_HINTS ${DIRECTORY_HINTS}")
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
            if(NOT TARGET OpenBLAS)
                add_library(OpenBLAS ${LINK_TYPE} IMPORTED)
                set_target_properties(OpenBLAS PROPERTIES
                        IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}"
                        INTERFACE_COMPILE_DEFINITIONS "OpenBLAS_AVAILABLE")
            endif()
        endif()
    endif()




    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS lib in system")
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
#            add_library(OpenBLAS INTERFACE IMPORTED)
#            target_link_libraries(OpenBLAS INTERFACE ${OpenBLAS_LIBRARIES})
#            target_include_directories(OpenBLAS INTERFACE ${OpenBLAS_INCLUDE_DIRS})
            add_library(OpenBLAS STATIC IMPORTED)
            set_target_properties(OpenBLAS PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}"
                    INTERFACE_COMPILE_DEFINITIONS "OpenBLAS_AVAILABLE"
                    )

        endif()
    endif()

    if(TARGET OpenBLAS)
        message(STATUS "OpenBLAS found: ${OpenBLAS_LIBRARIES}")
        message(STATUS "                ${OpenBLAS_INCLUDE_DIRS}")

        target_link_libraries(OpenBLAS INTERFACE gfortran pthread)
        #For convenience, define these variables
        add_library(blas   INTERFACE IMPORTED)
        add_library(lapack INTERFACE IMPORTED)
        target_link_libraries(blas   INTERFACE OpenBLAS)
        target_link_libraries(lapack INTERFACE OpenBLAS)
    else()
        message(STATUS "Could not find OpenBLAS")
    endif()
endfunction()