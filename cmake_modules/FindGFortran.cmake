
enable_language(Fortran)
if(GFORTRAN_LIB)
    return()
endif()

if(CMAKE_Fortran_COMPILER)
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX}
            OUTPUT_VARIABLE _libgfortran_path
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
endif()

if(EXISTS ${_libgfortran_path})
    get_filename_component(GFORTRAN_PATH ${_libgfortran_path} PATH)
    find_library(GFORTRAN_LIB gfortran PATHS ${GFORTRAN_PATH})
    message(STATUS "Found gfortran library:   ${GFORTRAN_LIB}")
else()
    # if libgfortran wasn't found at this point, the installation is probably broken
    # Let's try to find the library nonetheless.
    find_library(GFORTRAN_LIB gfortran)
endif()

if (NOT GFORTRAN_LIB)
    message(WARNING "gfortran library is required but could not be found")
endif ()

# also need libquadmath.a
find_library(QUADMATH_LIB quadmath PATHS ${GFORTRAN_PATH})
if (QUADMATH_LIB)
    list (APPEND GFORTRAN_LIB ${QUADMATH_LIB})
else ()
    message (WARNING "quadmath could not be found")
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
