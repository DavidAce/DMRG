
#enable_language(Fortran)

if(CMAKE_Fortran_COMPILER)
    if(BUILD_SHARED_LIBS)
        execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${CMAKE_SHARED_LIBRARY_SUFFIX}
                OUTPUT_VARIABLE GFORTRAN_LIB
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
    else()
        execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX}
                OUTPUT_VARIABLE GFORTRAN_LIB
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
    endif()
endif()

if(GFORTRAN_LIB)
    message(STATUS "Found gfortran library:   ${GFORTRAN_LIB}")
else()
    # if libgfortran wasn't found at this point, the installation is probably broken
    # Let's try to find the library nonetheless.
    find_library(GFORTRAN_LIB gfortran)
endif()

if (NOT GFORTRAN_LIB)
    message(FATAL_ERROR "gfortran library is required but could not be found")
endif ()

# also need libquadmath.a
get_filename_component(GFORTRAN_PATH ${GFORTRAN_LIB} PATH)
find_library(QUADMATH_LIB NAMES quadmath PATHS ${GFORTRAN_PATH})
if (NOT QUADMATH_LIB)
    message (FATAL_ERROR "quadmath could not be found")
endif ()

if (${CMAKE_HOST_APPLE})
    # also need -lgcc_ext.10.5
    find_library(GCC_EXT_LIB gcc_ext.10.5 PATHS ${GFORTRAN_PATH})
    if (GCC_EXT_LIB)
        list (APPEND GFORTRAN_LIB ${GCC_EXT_LIB})
    else ()
        message(STATUS "gcc_ext is required on MAC but could not be found")
    endif ()
endif()


add_library(gfortran INTERFACE)
target_link_libraries(gfortran INTERFACE ${GFORTRAN_LIB}  ${QUADMATH_LIB} )
