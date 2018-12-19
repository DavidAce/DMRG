
set(SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})

# Find from MKL
if(MKL_FOUND)
    # Try finding lapacke as module library
    message(STATUS "Attempting to find LAPACKE in loaded environment modules.")
    find_path(LAPACKE_INCLUDE_DIRS
            NAMES mkl_lapacke.h
            PATHS "${MKL_INCLUDE_DIR}"
            NO_DEFAULT_PATH
            )
    if(DEFINED LAPACKE_INCLUDE_DIRS)
        set(LAPACKE_FOUND TRUE)
        message(STATUS "Found LAPACKE in Intel MKL")
        target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
        return()
    endif()
endif()



#Check if it already exists from project libs directory
if(NOT LAPACKE_FOUND AND LAPACKE_FROM_OPENBLAS)
    set(LAPACKE_LIBRARY      ${BLAS_LIBRARIES})
    set(LAPACKE_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
    set(LAPECKE_FOUND TRUE)
    message(STATUS "Found LAPACKE in OpenBLAS")
    target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACKE_LIBRARY})
    target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
    return()
endif()


if (NOT LAPACKE_FOUND AND EXISTS "$ENV{LAPACKE_DIR}")
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
    if(LAPACKE_LIBRARY AND LAPACKE_INCLUDE_DIRS)
        set(LAPECKE_FOUND TRUE)
        message(STATUS "Found LAPACKE as a module")
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACKE_LIBRARY})
        target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
        return()
    endif()
endif()


# Find installed in system
if(NOT LAPACKE_FOUND)
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
    if(LAPACKE_LIBRARY AND LAPACKE_INCLUDE_DIRS)
        set(LAPECKE_FOUND TRUE)
        message(STATUS "Found LAPACKE in system")
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACKE_LIBRARY})
        target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
        return()
    endif()
endif()


if(LAPACKE_FOUND)
    message(STATUS "LAPACKE FOUND   : ${LAPACKE_LIBRARY}")
    message(STATUS "LAPACKE INCLUDE : ${LAPACKE_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(lapacke UNKNOWN IMPORTED)
    set(LAPACKE_LIBRARIES   ${LAPACKE_LIBRARY})

    set_target_properties(lapacke PROPERTIES
            IMPORTED_LOCATION               ${LAPACKE_LIBRARIES}
            INTERFACE_INCLUDE_DIRECTORIES   ${LAPACKE_INCLUDE_DIRS}
            )

    target_link_libraries(${PROJECT_NAME} PRIVATE lapacke)
    target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACKE_INCLUDE_DIRS})
    set(LAPACKE_INCLUDE_DIR ${LAPACKE_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Could not find package: LAPACKE")
endif()