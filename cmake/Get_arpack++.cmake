
function(find_arpackpp)

    find_package(Lapacke REQUIRED) # Lapacke needed by arpack++
    if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
        foreach (tgt BLAS::BLAS;LAPACK::LAPACK;LAPACKE::LAPACKE;gfortran::gfortran)
            if(NOT TARGET ${tgt})
                list(APPEND ARPACKPP_MISSING_TARGET ${tgt})
                mark_as_advanced(ARPACKPP_MISSING_TARGET)
            endif()
        endforeach()
        if(ARPACKPP_MISSING_TARGET)
            message(FATAL_ERROR "arpack++: dependencies missing [${ARPACKPP_MISSING_TARGET}]")
        endif()
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    set(arpack++_ROOT ${DMRG_DEPS_INSTALL_DIR})
    find_package(arpack++ ${REQUIRED})
    if(NOT arpack++_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "arpack++ will be installed into ${DMRG_DEPS_INSTALL_DIR}/arpack++ on first build.")
        include(cmake/InstallPackage.cmake)
        install_package(arpack++ "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(arpack++ REQUIRED)
    endif()
    if (NOT arpack++_FOUND OR NOT TARGET ARPACK::ARPACK++)
        message(FATAL_ERROR "arpack++ not found")
    endif()
endfunction()

find_arpackpp()


