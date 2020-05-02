find_package(Lapacke) # Lapacke needed by arpack++
if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    foreach (tgt blas::blas;lapack::lapack;lapacke::lapacke;gfortran::gfortran)
        if(NOT TARGET ${tgt})
            list(APPEND ARPACKPP_MISSING_TARGET ${tgt})
            mark_as_advanced(ARPACKPP_MISSING_TARGET)
        endif()
    endforeach()
    if(ARPACKPP_MISSING_TARGET)
        message(FATAL_ERROR "arpack++: dependencies missing [${ARPACKPP_MISSING_TARGET}]")
    endif()
endif()

if(NOT TARGET arpack::arpack++ AND DMRG_DOWNLOAD_METHOD STREQUAL "find")
    find_package(arpack++ REQUIRED)
    if (TARGET arpack::arpack++)
        message(STATUS "Found arpack++")
    endif()
endif()


if(NOT TARGET arpack::arpack++ AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    find_package(arpack++)
    if (TARGET arpack::arpack++)
        message(STATUS "Found arpack++")
    endif()
endif()

if(NOT TARGET arpack::arpack++ AND DMRG_DOWNLOAD_METHOD MATCHES "fetch")
    message(STATUS "arpack++ will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++ on first build.")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(arpack++ "${CMAKE_INSTALL_PREFIX}" "")
    find_package(arpack++ REQUIRED)
    if (TARGET arpack::arpack++)
        message(STATUS "Successfully installed arpack++")
    endif()
endif()



