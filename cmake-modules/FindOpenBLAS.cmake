function(openblas_message TYPE MSG)
    if(OPENBLAS_FIND_VERBOSE)
        message(${TYPE} ${MSG})
    endif()
endfunction()

function(find_OpenBLAS)

    # We can't use openblas from conda when compiling with Clang.
    # We seem to be able to use the one from apt though, so we add /usr to
    # the hints.
    # This means we are more likely to build it from source on Clang.
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()



    if(NOT TARGET openblas::openblas)
        openblas_message(STATUS "Looking for OpenBLAS config")
        find_package(OpenBLAS 0.3
                HINTS ${CMAKE_INSTALL_PREFIX}
                PATH_SUFFIXES
                openblas lib/x86_64-linux-gnu
                ${NO_DEFAULT_PATH} ${NO_CMAKE_PACKAGE_REGISTRY}
                CONFIG
                )


        # Some OpenBLAS cmake builds generate an OpenBLAS::OpenBLAS target. Ignore it.
        # Why? - It links "-lpthread" instead of "pthread" and some versions inject shared libraries
        # on static builds. Instead, make sure to import the correct static/shared libraries manually.
        if(TARGET OpenBLAS::OpenBLAS)
            get_target_property(OpenBLAS_LIBRARIES    OpenBLAS::OpenBLAS LOCATION)
            get_target_property(OpenBLAS_INCLUDE_DIRS OpenBLAS::OpenBLAS INTERFACE_INCLUDE_DIRECTORIES)
        endif()

        if(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            get_filename_component(OpenBLAS_LIBRARIES_EXT ${OpenBLAS_LIBRARIES} EXT)
            if(NOT "${OpenBLAS_LIBRARIES_EXT}" MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}" AND NOT BUILD_SHARED_LIBS)
                get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
                get_filename_component(OPENBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
                if(EXISTS  ${OPENBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX} )
                    openblas_message(STATUS "Replaced OpenBLAS extension from ${OpenBLAS_LIBRARIES_EXT} to ${CMAKE_STATIC_LIBRARY_SUFFIX}")
                    set(OpenBLAS_LIBRARIES ${OPENBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CMAKE_STATIC_LIBRARY_SUFFIX})
                else()
                    openblas_message(STATUS "Could not find OpenBLAS library with correct suffix (i.e. static/shared as requested): ${OpenBLAS_LIBRARIES}")
                endif()
            endif()
            if(NOT TARGET openblas::openblas)
                add_library(openblas::openblas ${LINK_TYPE} IMPORTED)
                set_target_properties(openblas::openblas PROPERTIES
                        IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                        INTERFACE_COMPILE_DEFINITIONS "OPENBLAS_AVAILABLE")
                target_include_directories(openblas::openblas SYSTEM INTERFACE ${OpenBLAS_INCLUDE_DIRS})

            endif()
        endif()
    endif()




    if(NOT TARGET openblas::openblas)
        openblas_message(STATUS "Looking for OpenBLAS in system")

        find_library(OpenBLAS_LIBRARIES
                NAMES openblas
                HINTS ${CMAKE_INSTALL_PREFIX}
                PATHS
                PATH_SUFFIXES
                lib openblas/lib OpenBLAS/lib openblas OpenBLAS lib/x86_64-linux-gnu
                ${NO_DEFAULT_PATH}
                ${NO_CMAKE_PACKAGE_REGISTRY}
                )
        find_path(OpenBLAS_INCLUDE_DIRS
                NAMES openblas_config.h
                HINTS ${CMAKE_INSTALL_PREFIX}
                PATH_SUFFIXES
                include openblas openblas/include OpenBLAS OpenBLAS/include blas/include include/x86_64-linux-gnu
                ${NO_DEFAULT_PATH}
                ${NO_CMAKE_PACKAGE_REGISTRY}
                )
        if (OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
            add_library(openblas::openblas ${LINK_TYPE} IMPORTED)
            set_target_properties(openblas::openblas PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                    INTERFACE_COMPILE_DEFINITIONS "OPENBLAS_AVAILABLE")
            target_include_directories(openblas::openblas SYSTEM INTERFACE ${OpenBLAS_INCLUDE_DIRS})
        endif()
    endif()

endfunction()


find_OpenBLAS()

if(TARGET openblas::openblas)
    get_target_property(OpenBLAS_LIBRARIES openblas::openblas LOCATION)
    get_target_property(OpenBLAS_INCLUDE_DIRS openblas::openblas INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    openblas_message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")

    set(OpenBLAS_TARGET openblas::openblas)
    target_link_libraries(openblas::openblas INTERFACE gfortran::gfortran pthread)

    # Fix for OpenBLAS 0.3.9, which otherwise includes <complex> inside of an extern "C" scope.
    target_compile_definitions(openblas::openblas INTERFACE LAPACK_COMPLEX_CUSTOM)
    target_compile_definitions(openblas::openblas INTERFACE lapack_complex_float=std::complex<float>)
    target_compile_definitions(openblas::openblas INTERFACE lapack_complex_double=std::complex<double>)

    #For convenience, define these targes
    add_library(blas::blas                  INTERFACE IMPORTED)
    add_library(lapack::lapack              INTERFACE IMPORTED)
    target_link_libraries(blas::blas        INTERFACE openblas::openblas)
    target_link_libraries(lapack::lapack    INTERFACE openblas::openblas)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
        OpenBLAS
        DEFAULT_MSG OpenBLAS_LIBRARIES OpenBLAS_INCLUDE_DIRS OpenBLAS_TARGET
        )