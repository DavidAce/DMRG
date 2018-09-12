
set(SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})

#Check if it already exists from project libs directory
if(LAPACKE_FROM_OPENBLAS)
    set(LAPACKE_LIBRARY      ${BLAS_LIBRARIES})
    set(LAPACKE_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
endif()


if (NOT LAPACKE_LIBRARY AND EXISTS "$ENV{LAPACKE_DIR}")
    # Try finding lapacke as module library
    message(STATUS "Attempting to find LAPACKE in loaded environment modules.")
    find_library(LAPACKE_LIBRARY
            NAMES libopenblas${SUFFIX}
            PATHS "$ENV{BLAS_DIR}/lib"
            NO_DEFAULT_PATH
            )

    find_path(LAPACKE_INCLUDE_DIRS
            NAMES lapacke.h
            PATHS "$ENV{BLAS_DIR}/include"
            NO_DEFAULT_PATH
            )
endif()

# Find installed in system
if(NOT LAPACKE_LIBRARY)
    message(STATUS "Attempting to find LAPACKE in system")
    find_library(LAPACKE_LIBRARY
            NAMES liblapacke${SUFFIX}
            PATHS "/usr/lib/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
    find_path(LAPACKE_INCLUDE_DIRS
            NAMES lapacke.h
            PATHS "/usr/include" "/usr/include/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
endif()


if(LAPACKE_LIBRARY)
    message(STATUS "LAPACKE FOUND   : ${LAPACKE_LIBRARY}")
    message(STATUS "LAPACKE INCLUDE : ${LAPACKE_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(lapacke UNKNOWN IMPORTED)
    set(LAPACKE_LIBRARIES   ${LAPACKE_LIBRARY})

    set_target_properties(lapacke PROPERTIES
            IMPORTED_LOCATION               "${LAPACKE_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORY     "${LAPACKE_INCLUDE_DIRS}"
            INTERFACE_LINK_FLAGS            ""
            INTERFACE_COMPILE_OPTIONS       ""
            )

    target_link_libraries(${PROJECT_NAME} PRIVATE lapacke)
    target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Could not find package: LAPACKE")
endif()