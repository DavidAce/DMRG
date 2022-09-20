function(find_gfortran)
    if(NOT BUILD_SHARED_LIBS)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_SUFFIX})
        set(LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    else()
        set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})

    endif()
    include(${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)

    if(CMAKE_Fortran_COMPILER)
        execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${LIB_SUFFIX}
                        OUTPUT_VARIABLE GFORTRAN_LIBRARY
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        )
        execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libquadmath${LIB_SUFFIX}
                        OUTPUT_VARIABLE QUADMATH_LIBRARY
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        )
    endif()
    if(IS_ABSOLUTE "${GFORTRAN_LIBRARY}" AND EXISTS "${GFORTRAN_LIBRARY}")
        set(GFORTRAN_LIBRARY "${GFORTRAN_LIBRARY}" CACHE PATH "Path to gfortran library")
    else()
        message(STATUS "GFORTRAN_LIBRARY does not exist: ${GFORTRAN_LIBRARY}")
        unset(GFORTRAN_LIBRARY)
        unset(GFORTRAN_LIBRARY CACHE)
        find_library(GFORTRAN_LIBRARY gfortran)
    endif()
    if(IS_ABSOLUTE "${QUADMATH_LIBRARY}" AND EXISTS "${QUADMATH_LIBRARY}")
        set(QUADMATH_LIBRARY "${QUADMATH_LIBRARY}" CACHE PATH "Path to quadmath library")
    else()
        unset(QUADMATH_LIBRARY)
        unset(QUADMATH_LIBRARY CACHE)
        find_library(QUADMATH_LIBRARY quadmath)
    endif()
endfunction()

find_gfortran()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(gfortran DEFAULT_MSG GFORTRAN_LIBRARY)

if(gfortran_FOUND)
    if(NOT TARGET gfortran::gfortran)
        add_library(gfortran::gfortran UNKNOWN IMPORTED)
        set_target_properties(gfortran::gfortran PROPERTIES IMPORTED_LOCATION "${GFORTRAN_LIBRARY}")
        if(QUADMATH_LIB)
            add_library(gfortran::quadmath UNKNOWN IMPORTED)
            set_target_properties(gfortran::quadmath PROPERTIES IMPORTED_LOCATION "${QUADMATH_LIBRARY}")
            target_link_libraries(gfortran::gfortran INTERFACE gfortran::quadmath)
        endif()
        add_library(gfortran ALIAS gfortran::gfortran)
    endif()
endif()