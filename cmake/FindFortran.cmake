function(find_gfortran)
    if(NOT TARGET gfortran::gfortran)
        #enable_language(Fortran)
        include(${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
        if(BUILD_SHARED_LIBS)
            set(GFORTRAN_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
        else()
            set(GFORTRAN_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
        endif()

        if(CMAKE_Fortran_COMPILER)
            execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${GFORTRAN_LIB_SUFFIX}
                    OUTPUT_VARIABLE GFORTRAN_LIB
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
            execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libquadmath${GFORTRAN_LIB_SUFFIX}
                    OUTPUT_VARIABLE QUADMATH_LIB
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
        endif()


        if(GFORTRAN_LIB)
            message(STATUS "Found gfortran library:   ${GFORTRAN_LIB}")
        else()
            # if libgfortran wasn't found at this point, the installation is probably broken
            # Let's try to find the library nonetheless.
            find_library(GFORTRAN_LIB libgfortran${GFORTRAN_LIB_SUFFIX})
        endif()
        if (NOT GFORTRAN_LIB)
            message(FATAL_ERROR "gfortran library is required but could not be found")
        endif ()

        if(QUADMATH_LIB)
            message(STATUS "Found quadmath library: ${QUADMATH_LIB}")
        else()
            # if libgfortran wasn't found at this point, the installation is probably broken
            # Let's try to find the library nonetheless.
            find_library(QUADMATH_LIB libquadmath${GFORTRAN_LIB_SUFFIX})
        endif()
        if (NOT QUADMATH_LIB)
            message(FATAL_ERROR "quadmath library is required but could not be found")
        endif ()


        if (${CMAKE_HOST_APPLE})
            # also need -lgcc_ext.10.5
            find_library(GCC_EXT_LIB gcc_ext.10.5 PATHS ${GFORTRAN_PATH})
            if (GCC_EXT_LIB)
                list (APPEND GFORTRAN_LIB ${GCC_EXT_LIB})
            else ()
                message(WARNING "gcc_ext is required on MAC but could not be found")
            endif ()
        endif()

        add_library(gfortran_lib ${LINK_TYPE} IMPORTED)
        add_library(quadmath_lib ${LINK_TYPE} IMPORTED)
        set_target_properties(quadmath_lib PROPERTIES IMPORTED_LOCATION "${QUADMATH_LIB}")
        set_target_properties(gfortran_lib PROPERTIES IMPORTED_LOCATION "${GFORTRAN_LIB}" INTERFACE_LINK_LIBRARIES quadmath_lib)
        add_library(gfortran::gfortran IMPORTED INTERFACE)
        target_link_libraries(gfortran::gfortran INTERFACE quadmath_lib gfortran_lib)
        mark_as_advanced(GFORTRAN_LIB_SUFFIX)
        mark_as_advanced(GFORTRAN_LIB)
        mark_as_advanced(QUADMATH_LIB)
    endif()
endfunction()


find_gfortran()
if(TARGET gfortran::gfortran)
    get_target_property(GFORTRAN_LIB gfortran_lib LOCATION)
    get_target_property(QUADMATH_LIB quadmath_lib LOCATION)
    set(GFORTRAN_TARGET gfortran::gfortran)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Fortran
        DEFAULT_MSG GFORTRAN_LIB QUADMATH_LIB GFORTRAN_TARGET
        )
