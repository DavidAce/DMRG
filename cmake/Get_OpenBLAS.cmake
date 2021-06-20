
# If the INTEL MKL library has already been loaded, skip the rest.
if(TARGET mkl::mkl)
    return()
endif()

function(find_openblas)
    if(DMRG_PACKAGE_MANAGER MATCHES "find")
        set(REQUIRED REQUIRED)
    endif()
    find_package(OpenBLAS
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            COMPONENTS ${OpenBLAS_LINKAGE} ${OpenBLAS_THREADING}
            ${REQUIRED})

    if(NOT OpenBLAS_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "OpenBLAS will be installed into ${DMRG_DEPS_INSTALL_DIR}/OpenBLAS")
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DTARGET:STRING=${OPENBLAS_MARCH})
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DUSE_THREAD:BOOL=1)
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DBUILD_RELAPACK:BOOL=OFF)
        include(cmake/InstallPackage.cmake)
        install_package(OpenBLAS "${DMRG_DEPS_INSTALL_DIR}" "${OPENBLAS_CMAKE_OPTIONS}")
        find_package(OpenBLAS
                HINTS ${DMRG_DEPS_INSTALL_DIR}
                COMPONENTS ${OpenBLAS_LINKAGE} ${OpenBLAS_THREADING}
                REQUIRED)
    endif()
    if(OpenBLAS_FOUND)
        message(STATUS "Found OpenBLAS: ${OpenBLAS_FOUND}")
    endif()
    if(NOT TARGET OpenBLAS::OpenBLAS)
        message(FATAL_ERROR "Target undefined: OpenBLAS::OpenBLAS")
    endif()
    #For convenience, define these targes
    if(NOT TARGET BLAS::BLAS)
        add_library(BLAS::BLAS                  INTERFACE IMPORTED)
        target_link_libraries(BLAS::BLAS        INTERFACE OpenBLAS::OpenBLAS)
    endif()
    if(NOT TARGET LAPACK::LAPACK)
        add_library(LAPACK::LAPACK              INTERFACE IMPORTED)
        target_link_libraries(LAPACK::LAPACK    INTERFACE OpenBLAS::OpenBLAS)
    endif()
endfunction()

find_openblas()


