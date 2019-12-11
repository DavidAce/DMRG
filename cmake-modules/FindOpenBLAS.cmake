function(find_OpenBLAS)
    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS config")
        find_package(OpenBLAS 0.3
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{EBROOTBLAS}
                    $ENV{OpenBLAS_DIR}      ${OpenBLAS_DIR}
                    $ENV{openblas_DIR}      ${openblas_DIR}
                    $ENV{BLAS_DIR}          ${BLAS_DIR}
                    $ENV{blas_DIR}          ${blas_DIR}
                    $ENV{LAPACKE_DIR}       ${LAPACKE_DIR}
                    $ENV{lapacke_DIR}       ${lapacke_DIR}
                PATH_SUFFIXES
                    lib OpenBLAS/lib OpenBLAS openblas
                )

        if(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
            get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
            if(EXISTS  ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
                set(OpenBLAS_LIBRARIES ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
                add_library(OpenBLAS INTERFACE IMPORTED)
                target_link_libraries(OpenBLAS INTERFACE ${OpenBLAS_LIBRARIES})
                target_include_directories(OpenBLAS INTERFACE ${OpenBLAS_INCLUDE_DIRS})
            else()
                message(STATUS "Could not find OpenBLAS library with correct suffix (i.e. static/shared as requested): ${OpenBLAS_LIBRARIES}")
            endif()
        endif()
    endif()




    if(NOT TARGET OpenBLAS)
        message(STATUS "Searching for OpenBLAS lib in system")
        find_library(OpenBLAS_LIBRARIES
                NAMES libopenblas${CUSTOM_SUFFIX}
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{EBROOTBLAS}
                    $ENV{OpenBLAS_DIR}      ${OpenBLAS_DIR}
                    $ENV{openblas_DIR}      ${openblas_DIR}
                    $ENV{BLAS_DIR}          ${BLAS_DIR}
                    $ENV{blas_DIR}          ${blas_DIR}
                    $ENV{LAPACKE_DIR}       ${LAPACKE_DIR}
                    $ENV{lapacke_DIR}       ${lapacke_DIR}
                PATH_SUFFIXES
                    lib openblas/lib OpenBLAS/lib openblas OpenBLAS

                )
        find_path(OpenBLAS_INCLUDE_DIRS
                NAMES openblas_config.h
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{EBROOTBLAS}
                    $ENV{OpenBLAS_DIR}    ${OpenBLAS_DIR}
                    $ENV{openblas_DIR}    ${openblas_DIR}
                    $ENV{BLAS_DIR}    ${BLAS_DIR}
                    $ENV{blas_DIR}    ${blas_DIR}
                    $ENV{LAPACKE_DIR} ${LAPACKE_DIR}
                PATH_SUFFIXES
                    include openblas openblas/include OpenBLAS OpenBLAS/include blas/include
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