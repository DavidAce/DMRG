function(lapacke_message TYPE MSG)
    if(LAPACKE_FIND_VERBOSE)
        message(${TYPE} ${MSG})
    endif()
endfunction()



function(find_Lapacke)
    include(cmake-modules/CheckLapackeCompiles.cmake)
    if (NOT TARGET lapacke::lapacke)
        # Find from MKL
        if(TARGET mkl::mkl)
            # Try finding lapacke in MKL library
            lapacke_message(STATUS "Looking for Lapacke in Intel MKL")
            if(MKL_INCLUDE_DIR)
                CheckLapackeCompiles("mkl::mkl" ""  ""  ""  "")
            endif()

            if(LAPACKE_COMPILES)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                get_target_property(LAPACKE_INCLUDE_DIR mkl::mkl INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
                target_link_libraries(lapacke::lapacke INTERFACE mkl::mkl)
                target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                lapacke_message(STATUS "Looking for Lapacke in Intel MKL - found")
            else()
                message(STATUS "Looking for Lapacke in Intel MKL - not found")
                message(WARNING "Intel MKL is expected to bundle Lapacke")
            endif()
        endif()
    endif()

#
    if (NOT TARGET lapacke::lapacke)
        if (TARGET openblas::openblas)
            lapacke_message(STATUS "Looking for Lapacke in OpenBLAS")
            CheckLapackeCompiles("openblas::openblas" "" "" "" "")
            if(LAPACKE_COMPILES)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                get_target_property(LAPACKE_INCLUDE_DIR openblas::openblas INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
                target_link_libraries(lapacke::lapacke INTERFACE openblas::openblas)
                target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                lapacke_message(STATUS "Looking for Lapacke in OpenBLAS - found")
            else()
                message(STATUS "Looking for Lapacke in OpenBLAS - not found")
                message(WARNING "OpenBLAS is expected to bundle Lapacke")
            endif()
        endif()
    endif()



    if (NOT TARGET lapacke::lapacke)
        if(TARGET lapack::lapack)
            lapacke_message(STATUS "Looking for Lapacke in system")
            find_library(LAPACKE_LIBRARY
                    NAMES lapacke
                    PATH_SUFFIXES
                    OpenBLAS openblas openblas/lib OpenBLAS/lib lapack lapack/lib blas blas/lib
                    )
            find_path(LAPACKE_INCLUDE_DIR
                    NAMES lapacke.h
                    PATH_SUFFIXES
                        OpenBLAS openblas openblas/include OpenBLAS/include lapack)

            if(NOT LAPACKE_LIBRARY)
                set(LAPACKE_LIBRARY "")
            endif()
            if(NOT LAPACKE_INCLUDE_DIR)
                set(LAPACKE_INCLUDE_DIR "")
            endif()
            CheckLapackeCompiles("lapack::lapack" "${LAPACKE_LIBRARY}" "${LAPACKE_INCLUDE_DIR}" "" "")
            if(LAPACKE_COMPILES)
                if(LAPACKE_LIBRARY)
                    add_library(lapacke::lapacke ${LINK_TYPE} IMPORTED)
                    set_target_properties(lapacke::lapacke PROPERTIES IMPORTED_LOCATION "${LAPACKE_LIBRARY}")
                else()
                    add_library(lapacke::lapacke INTERFACE IMPORTED)
                endif()
                target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
            endif()

            if(TARGET lapacke::lapacke)
                target_link_libraries(lapacke::lapacke INTERFACE blas::blas lapack::lapack)
                lapacke_message(STATUS "Looking for Lapacke in system - found: ${LAPACKE_LIBRARY}")
            else()
                message(STATUS "Looking for Lapacke in system - not found")
                message(STATUS "Tried LAPACKE_LIBRARY     : ${LAPACKE_LIBRARY} ")
                message(STATUS "Tried LAPACKE_INCLUDE_DIR : ${LAPACKE_INCLUDE_DIR} ")
                message(WARNING "No version of lapacke could be found")
            endif()
        endif()
    endif()

endfunction()

find_Lapacke()
if(TARGET lapacke::lapacke)
    get_target_property(LAPACKE_INCLUDE_DIR lapacke::lapacke INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    set(LAPACKE_TARGET lapacke::lapacke)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapacke
        DEFAULT_MSG LAPACKE_INCLUDE_DIR LAPACKE_TARGET
)