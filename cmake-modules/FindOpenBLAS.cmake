function(find_OpenBLAS)

    # We can't use openblas from conda when compiling with Clang.
    # We seem to be able to use the one from apt though, so we add /usr to
    # the hints.
    # This means we are more likely to build it from source on Clang.
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        set(OPENBLAS_HINTS /usr ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS})
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    else()
        set(OPENBLAS_HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS})
    endif()



    if(NOT TARGET openblas::openblas)
        if(OpenBLAS_FIND_VERBOSE)
            message(STATUS "Searching for OpenBLAS config")
        endif()
        find_package(OpenBLAS 0.3
                HINTS ${OPENBLAS_HINTS}
                PATHS
                $ENV{EBROOTBLAS}
                $ENV{BLAS_DIR}
                $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                lib OpenBLAS/lib OpenBLAS openblas lib/x86_64-linux-gnu
                ${NO_DEFAULT_PATH}
                QUIET
                )

        if(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            get_filename_component(OpenBLAS_LIBRARIES_EXT ${OpenBLAS_LIBRARIES} EXT)
            if(NOT "${OpenBLAS_LIBRARIES_EXT}" MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}" AND NOT BUILD_SHARED_LIBS)
                get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
                get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
                if(EXISTS  ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX} )
                    if(OpenBLAS_FIND_VERBOSE)
                        message(STATUS "Replaced OpenBLAS extension from ${OpenBLAS_LIBRARIES_EXT} to ${CMAKE_STATIC_LIBRARY_SUFFIX}")
                    endif()
                    set(OpenBLAS_LIBRARIES ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX})
                else()
                    if(OpenBLAS_FIND_VERBOSE)
                        message(STATUS "Could not find OpenBLAS library with correct suffix (i.e. static/shared as requested): ${OpenBLAS_LIBRARIES}")
                    endif()
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
        if(OpenBLAS_FIND_VERBOSE)
            message(STATUS "Searching for OpenBLAS in system")
        endif()
        find_library(OpenBLAS_LIBRARIES
                NAMES openblas
                HINTS ${OPENBLAS_HINTS}
                PATHS
                $ENV{EBROOTBLAS}
                $ENV{BLAS_DIR}
                $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                lib openblas/lib OpenBLAS/lib openblas OpenBLAS lib/x86_64-linux-gnu
                ${NO_DEFAULT_PATH}
                )
        find_path(OpenBLAS_INCLUDE_DIRS
                NAMES openblas_config.h
                HINTS ${OPENBLAS_HINTS}
                PATHS
                $ENV{EBROOTBLAS}
                $ENV{BLAS_DIR}
                $ENV{BLAS_ROOT}
                PATH_SUFFIXES
                include openblas openblas/include OpenBLAS OpenBLAS/include blas/include include/x86_64-linux-gnu
                ${NO_DEFAULT_PATH}
                )
        if (OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            add_library(openblas::openblas STATIC IMPORTED)
            set_target_properties(openblas::openblas PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                    INTERFACE_COMPILE_DEFINITIONS "OpenBLAS_AVAILABLE")
            target_include_directories(openblas::openblas SYSTEM INTERFACE ${OpenBLAS_INCLUDE_DIRS})
        endif()
    endif()

endfunction()


find_OpenBLAS()
if(TARGET openblas::openblas)
    message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
    get_target_property(OpenBLAS_LIBRARIES openblas::openblas LOCATION)
    get_target_property(OpenBLAS_INCLUDE_DIRS openblas::openblas INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    set(OpenBLAS_TARGET openblas::openbla)
    set(OpenBLAS_FOUND  TRUE)
    target_link_libraries(openblas::openblas INTERFACE gfortran::gfortran pthread)
    #For convenience, define these targes
    add_library(blas::blas                  INTERFACE IMPORTED)
    add_library(lapack::lapack              INTERFACE IMPORTED)
    target_link_libraries(blas::blas        INTERFACE openblas::openblas)
    target_link_libraries(lapack::lapack    INTERFACE openblas::openblas)
else()
    set(OpenBLAS_FOUND FALSE)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
        OpenBLAS
        FOUND_VAR OpenBLAS_FOUND
        REQUIRED_VARS OpenBLAS_LIBRARIES OpenBLAS_INCLUDE_DIRS OpenBLAS_TARGET
        HANDLE_COMPONENTS
        FAIL_MESSAGE "Failed to find OpenBLAS"
        )