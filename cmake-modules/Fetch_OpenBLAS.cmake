
# If the INTEL MKL library has already been loaded, skip the rest.
if(TARGET mkl::mkl)
    return()
endif()


if(NOT TARGET openblas::openblas AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    set(OpenBLAS_FIND_VERBOSE ON)
    find_package(OpenBLAS)
endif()

if(NOT TARGET openblas::openblas AND DMRG_DOWNLOAD_METHOD MATCHES "fetch|native")
    message(STATUS "OpenBLAS will be installed into ${CMAKE_INSTALL_PREFIX}/OpenBLAS")
    set(OPENBLAS_MULTITHREADED 1 )
    if(TARGET openmp::openmp AND BUILD_SHARED_LIBS)
        set(OPENBLAS_ENABLE_OPENMP 1)
    else()
        set(OPENBLAS_ENABLE_OPENMP 0) # OpenMP doesnt work on static builds with OpenBLAS
    endif()
    include(cmake-modules/getExpandedTarget.cmake)
    #expand_target_libs(gfortran::gfortran GFORTRAN_LIBS)
    #list(FILTER GFORTRAN_LIBS INCLUDE REGEX "gfortran|quadmath")
    #list(JOIN GFORTRAN_LIBS " " GFORTRAN_LIBS_SPACE_SEP_STRING)
    #list(GET  GFORTRAN_LIBS 0 GFORTRAN_LIB)
    #get_filename_component(GFORTRAN_PATH ${GFORTRAN_LIB} DIRECTORY)

    #set(LDFLAGS "-L${GFORTRAN_PATH} ${GFORTRAN_LIBS_SPACE_SEP_STRING}")
    #set(FFLAGS  "-O3 -Wno-maybe-uninitialized -Wno-conversion -Wno-unused-but-set-variable -Wno-unused-variable")
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
