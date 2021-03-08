function(lapacke_message TYPE MSG)
    if(NOT Lapacke_FIND_QUIETLY)
        message(${TYPE} ${MSG})
    endif()
endfunction()



function(find_Lapacke)
    include(cmake-modules/CheckLapackeCompiles.cmake)
    if (NOT TARGET LAPACKE::LAPACKE)
        # Find from MKL
        if(TARGET mkl::mkl)
            # Try finding lapacke in MKL library
            lapacke_message(STATUS "Looking for Lapacke in Intel MKL")
            if(MKL_INCLUDE_DIR)
                check_lapacke_compiles("mkl::mkl" ""  ""  ""  "")
            endif()

            if(LAPACKE_COMPILES)
                add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
                target_link_libraries(LAPACKE::LAPACKE INTERFACE mkl::mkl)
                get_target_property(MKL_INCLUDE_DIR mkl::mkl INTERFACE_INCLUDE_DIRECTORIES)
                if(MKL_INCLUDE_DIR)
                    target_include_directories(LAPACKE::LAPACKE SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
                    set(LAPACKE_INCLUDE_DIR ${MKL_INCLUDE_DIR})
                endif()
                mark_as_advanced(MKL_INCLUDE_DIR)
                lapacke_message(STATUS "Looking for Lapacke in Intel MKL - found: ${LAPACKE_INCLUDE_DIR}")
            else()
                message(STATUS "Looking for Lapacke in Intel MKL - not found")
                message(WARNING "Intel MKL is expected to bundle Lapacke")
            endif()
        endif()
    endif()

#
    if (NOT TARGET LAPACKE::LAPACKE)
        if (TARGET openblas::openblas)
            lapacke_message(STATUS "Looking for Lapacke in OpenBLAS")
            check_lapacke_compiles("openblas::openblas" "" "" "" "")
            if(LAPACKE_COMPILES)
                add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
                get_target_property(LAPACKE_INCLUDE_DIR openblas::openblas INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
                target_link_libraries(LAPACKE::LAPACKE INTERFACE openblas::openblas)
                target_include_directories(LAPACKE::LAPACKE SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                lapacke_message(STATUS "Looking for Lapacke in OpenBLAS - found")
            else()
                message(STATUS "Looking for Lapacke in OpenBLAS - not found")
                message(WARNING "OpenBLAS is expected to bundle Lapacke")
            endif()
        endif()
    endif()



    if (NOT TARGET LAPACKE::LAPACKE)
        if(TARGET LAPACK::LAPACK)
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
            check_lapacke_compiles("LAPACK::LAPACK" "${LAPACKE_LIBRARY}" "${LAPACKE_INCLUDE_DIR}" "" "")
            if(LAPACKE_COMPILES)
                if(LAPACKE_LIBRARY)
                    add_library(LAPACKE::LAPACKE ${LINK_TYPE} IMPORTED)
                    set_target_properties(LAPACKE::LAPACKE PROPERTIES IMPORTED_LOCATION "${LAPACKE_LIBRARY}")
                else()
                    add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
                endif()
                target_include_directories(LAPACKE::LAPACKE SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
            endif()

            if(TARGET LAPACKE::LAPACKE)
                target_link_libraries(LAPACKE::LAPACKE INTERFACE BLAS::BLAS LAPACK::LAPACK)
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
if(TARGET LAPACKE::LAPACKE)
    get_target_property(LAPACKE_INCLUDE_DIR LAPACKE::LAPACKE INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    set(LAPACKE_TARGET LAPACKE::LAPACKE)
else()
    unset(LAPACKE_COMPILES)
    unset(LAPACKE_COMPILES CACHE)
    unset(LAPACKE_COMPILES PARENT_SCOPE)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapacke
        DEFAULT_MSG LAPACKE_INCLUDE_DIR LAPACKE_TARGET
)