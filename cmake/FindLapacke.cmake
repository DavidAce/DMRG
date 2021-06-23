function(lapacke_message TYPE MSG)
    if(NOT Lapacke_FIND_QUIETLY)
        message(${TYPE} ${MSG})
    endif()
endfunction()



function(find_Lapacke)
    include(cmake/CheckLapackeCompiles.cmake)
    set(lapacke_tgt_candidates BLAS::BLAS LAPACK::LAPACK mkl::mkl OpenBLAS::OpenBLAS)

    foreach(tgt ${lapacke_tgt_candidates})
        if(TARGET ${tgt})
            lapacke_message(STATUS "Looking for Lapacke in ${tgt}")
            check_lapacke_compiles("${tgt}" "" "" "" "")
            if(LAPACKE_COMPILES)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                target_link_libraries(lapacke::lapacke INTERFACE ${tgt})
                lapacke_message(STATUS "Looking for Lapacke in ${tgt} - found")
                break()
            endif()
        endif()
    endforeach()
endfunction()

find_Lapacke()
if(TARGET lapacke::lapacke)
    set(LAPACKE_TARGET lapacke::lapacke)
else()
    unset(LAPACKE_COMPILES)
    unset(LAPACKE_COMPILES CACHE)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapacke
        DEFAULT_MSG LAPACKE_TARGET LAPACKE_COMPILES
)