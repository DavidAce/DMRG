
# If the INTEL MKL library has already been loaded, skip the rest.
if(TARGET mkl::mkl)
    return()
endif()


if(NOT TARGET openblas::openblas AND DMRG_PACKAGE_MANAGER MATCHES "find")
    set(OpenBLAS_FIND_VERBOSE ON)
    find_package(OpenBLAS)
endif()

if(NOT TARGET openblas::openblas AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
    message(STATUS "OpenBLAS will be installed into ${CMAKE_INSTALL_PREFIX}/OpenBLAS")
    set(OPENBLAS_MULTITHREADED 1 )
    if(TARGET openmp::openmp AND BUILD_SHARED_LIBS)
        set(OPENBLAS_ENABLE_OPENMP 1)
    else()
        set(OPENBLAS_ENABLE_OPENMP 0) # OpenMP doesnt work on static builds with OpenBLAS
    endif()
    include(cmake-modules/getExpandedTarget.cmake)
    list(APPEND OPENBLAS_CMAKE_OPTIONS  -DTARGET:STRING=${OPENBLAS_MARCH})
    list(APPEND OPENBLAS_CMAKE_OPTIONS  -DUSE_THREAD:BOOL=${OPENBLAS_MULTITHREADED})
    list(APPEND OPENBLAS_CMAKE_OPTIONS  -DBUILD_RELAPACK:BOOL=OFF)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(OpenBLAS "${CMAKE_INSTALL_PREFIX}" "${OPENBLAS_CMAKE_OPTIONS}")
    find_package(OpenBLAS REQUIRED)
    if(TARGET openblas::openblas)
        message(STATUS "Successfully installed OpenBLAS")
    endif()
endif()
