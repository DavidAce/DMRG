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
    endif()
    if(IS_ABSOLUTE "${GFORTRAN_LIBRARY}" AND EXISTS "${GFORTRAN_LIBRARY}")
        set(GFORTRAN_LIBRARY "${GFORTRAN_LIBRARY}" CACHE PATH "Path to gfortran library")
    else()
        message(STATUS "GFORTRAN_LIBRARY does not exist: ${GFORTRAN_LIBRARY}")
        unset(GFORTRAN_LIBRARY)
        unset(GFORTRAN_LIBRARY CACHE)
        find_library(GFORTRAN_LIBRARY gfortran)
    endif()
endfunction()

find_gfortran()

if(quadmath IN_LIST gfortran_FIND_COMPONENTS)
    find_package(quadmath ${gfortran_FIND_REQUIRED_quadmath})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(gfortran DEFAULT_MSG GFORTRAN_LIBRARY)

if(gfortran_FOUND AND GFORTRAN_LIBRARY)
    if(NOT TARGET gfortran::gfortran)
        add_library(gfortran::gfortran UNKNOWN IMPORTED)
        message(DEBUG "Defined target gfortran::gfortran for library: ${GFORTRAN_LIBRARY}")
    endif()
    set_target_properties(gfortran::gfortran PROPERTIES IMPORTED_LOCATION "${GFORTRAN_LIBRARY}")
    if(quadmath IN_LIST gfortran_FIND_COMPONENTS AND TARGET quadmath::quadmath)
        target_link_libraries(gfortran::gfortran INTERFACE quadmath::quadmath)
    endif()
endif()