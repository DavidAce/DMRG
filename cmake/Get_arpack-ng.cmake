function(find_arpack_ng)

    if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
        foreach (tgt BLAS::BLAS;LAPACK::LAPACK;gfortran::gfortran)
            if(NOT TARGET ${tgt})
                list(APPEND arpack_ng_MISSING_TARGET ${tgt})
                mark_as_advanced(arpack_ng_MISSING_TARGET)
            endif()
        endforeach()
        if(arpack_ng_MISSING_TARGET)
            message(FATAL_ERROR "arpack-ng: dependencies missing [${arpack_ng_MISSING_TARGET}]")
        endif()
    endif()

    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "cmake")
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()
    find_package(arpack-ng 3.8.0
            HINTS ${arpack-ng_ROOT} ${DMRG_DEPS_INSTALL_DIR}
            ${NO_DEFAULT_PATH}
            ${REQUIRED})
    if(NOT arpack-ng_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "arpack-ng will be installed into ${DMRG_DEPS_INSTALL_DIR}/arpack-ng")
        include(cmake/InstallPackage.cmake)
        install_package(arpack-ng "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(arpack-ng 3.8.0
                HINTS ${arpack-ng_ROOT} ${DMRG_DEPS_INSTALL_DIR}
                NO_DEFAULT_PATH
                REQUIRED)
    endif()

    if(arpack-ng_FOUND AND TARGET ARPACK::ARPACK)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    else()
        message(FATAL_ERROR "arpack-ng not found")
    endif()
endfunction()

find_arpack_ng()