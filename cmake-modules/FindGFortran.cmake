
enable_language(Fortran)
add_library(gfortran STATIC IMPORTED)
#if(GFORTRAN_LIB AND QUADMATH_LIB)
#    return()
#endif()

if(CMAKE_Fortran_COMPILER)
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX}
            OUTPUT_VARIABLE _libgfortran_path
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
endif()

if(EXISTS ${_libgfortran_path})
    get_filename_component(GFORTRAN_PATH ${_libgfortran_path} PATH)
    find_library(GFORTRAN_LIB NAMES libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX} PATHS ${GFORTRAN_PATH})
    message(STATUS "Found gfortran library:   ${GFORTRAN_LIB}")
else()
    # if libgfortran wasn't found at this point, the installation is probably broken
    # Let's try to find the library nonetheless.
    find_library(GFORTRAN_LIB libgfortran.so)
endif()

if (NOT GFORTRAN_LIB)
    message(FATAL_ERROR "gfortran library is required but could not be found")
endif ()

# also need libquadmath.a
find_library(QUADMATH_LIB NAMES libquadmath.a PATHS ${GFORTRAN_PATH})
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

set_target_properties(gfortran PROPERTIES
        IMPORTED_LOCATION        "${GFORTRAN_LIB}"
        INTERFACE_LINK_LIBRARIES "${QUADMATH_LIB}")
#target_link_libraries(${PROJECT_NAME} PRIVATE gfortran)
