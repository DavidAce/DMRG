
# If the INTEL MKL library has already been loaded, skip the rest.
if(TARGET mkl::mkl)
    return()
endif()

function(find_openblas)
    find_package(Fortran REQUIRED)
    find_package(Threads REQUIRED)

    find_package(OpenBLAS 0.3.8 ${N5} ${N6} ${N7} ${N8} ${REQUIRED}) # Flags ignore system packages. See cmake/SetupPaths.cmake

    if(NOT OpenBLAS_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "OpenBLAS will be installed into ${DMRG_DEPS_INSTALL_DIR}/OpenBLAS")
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DTARGET:STRING=${OPENBLAS_MARCH})
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DUSE_THREAD:BOOL=1)
        list(APPEND OPENBLAS_CMAKE_OPTIONS  -DBUILD_RELAPACK:BOOL=OFF)
        include(cmake/InstallPackage.cmake)
        install_package(OpenBLAS "${DMRG_DEPS_INSTALL_DIR}" "${OPENBLAS_CMAKE_OPTIONS}")
        find_package(OpenBLAS 0.3.8 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()
    if(OpenBLAS_FOUND)
        message(STATUS "Found OpenBLAS: ${OpenBLAS_FOUND}")
    endif()
    if(TARGET OpenBLAS::OpenBLAS)
        get_target_property(OpenBLAS_LIBRARIES OpenBLAS::OpenBLAS LOCATION)
        get_target_property(OpenBLAS_INCLUDE_DIRS OpenBLAS::OpenBLAS INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
        target_link_libraries(OpenBLAS::OpenBLAS INTERFACE gfortran::gfortran Threads::Threads)
        # Fix for OpenBLAS 0.3.9, which otherwise includes <complex> inside of an extern "C" scope.
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        #    target_compile_definitions(openblas::openblas INTERFACE LAPACK_COMPLEX_CUSTOM)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
    else()
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


