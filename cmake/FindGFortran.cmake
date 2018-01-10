# This code to find libgfortran is taken from
#   https://github.com/m-reuter/arpackpp/blob/master/CMakeLists.txt



# Find libgfortran (static preferred)
IF(NOT WIN32)
    # Query gfortran to get the libgfortran.so path
    FIND_PROGRAM(_GFORTRAN_EXECUTABLE NAMES gfortran f95 gfortran-mp-4.9)
    IF(_GFORTRAN_EXECUTABLE)
        EXECUTE_PROCESS(COMMAND ${_GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.a
                OUTPUT_VARIABLE _libgfortran_path
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
    ENDIF()
    unset(_GFORTRAN_EXECUTABLE CACHE)

    IF(EXISTS ${_libgfortran_path})
        get_filename_component(GFORTRAN_PATH ${_libgfortran_path} PATH)
        find_library(GFORTRAN_LIB gfortran PATHS ${GFORTRAN_PATH})
        message("FOUND GFORTRAN :   ${GFORTRAN_LIB}")
    ELSE()
        # if libgfortran wasn't found at this point, the installation is probably broken
        # Let's try to find the library nonetheless.
        FIND_LIBRARY(GFORTRAN_LIB gfortran)
    ENDIF()
    IF (NOT GFORTRAN_LIB)
        MESSAGE(STATUS "gfortran is required but could not be found")
        SET(ABORT_CONFIG TRUE)
    ENDIF (NOT GFORTRAN_LIB)

    # also need libquadmath.a
    find_library(QUADMATH_LIB quadmath PATHS ${GFORTRAN_PATH})
    if (QUADMATH_LIB)
        list (APPEND GFORTRAN_LIB ${QUADMATH_LIB})
    else ()
        message (STATUS "quadmath could not be found")
        set (ABORT_CONFIG TRUE)
    endif ()

    IF (APPLE)
        # also need -lgcc_ext.10.5
        find_library(GCC_EXT_LIB gcc_ext.10.5 PATHS ${GFORTRAN_PATH})
        if (GCC_EXT_LIB)
            list (APPEND GFORTRAN_LIB ${GCC_EXT_LIB})
        else ()
            MESSAGE(STATUS "gcc_ext is required on MAC but could not be found")
            SET(ABORT_CONFIG TRUE)
        endif ()
    ENDIF(APPLE)

ENDIF(NOT WIN32)