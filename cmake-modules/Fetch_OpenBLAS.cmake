
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

# It seems that libopenblas from "apt" includes lapack in newer versions of Ubuntu. Trusty LTS 14.04 does not
# libopenblas bundled with lapack. Therefore we test if both lapack and blas are present to distinghuish these cases.
# Otherwise, arpack-ng will complain later about undefined references.


message(STATUS "SEARCHING FOR OpenBLAS IN SYSTEM...")
set(BLA_VENDOR OpenBLAS)
set(BLAS_VERBOSE OFF)
set(BLA_STATIC ${STATIC_BUILD})


if (EXISTS "$ENV{BLAS_DIR}")
    # Try finding openblas as module library
    message(STATUS "Attempting to find BLAS and LAPACK from environment modules.")
    find_library(BLAS_openblas_LIBRARY
            NAMES libopenblas${CUSTOM_SUFFIX}
            PATHS "$ENV{BLAS_DIR}/lib"
            NO_DEFAULT_PATH
            )

    find_path(BLAS_INCLUDE_DIRS
            NAMES openblas_config.h
            PATHS "$ENV{BLAS_DIR}/include"
            NO_DEFAULT_PATH
            )
    if (BLAS_openblas_LIBRARY)
        set(OpenBLAS_MULTITHREADED 1 )
    endif()
else()
    find_library(BLAS_openblas_LIBRARY
            NAMES libopenblas${CUSTOM_SUFFIX}
            PATHS "/usr/lib/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
    find_path(BLAS_INCLUDE_DIRS
            NAMES cblas.h
            PATHS "/usr/include" "/usr/include/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
    if (BLAS_openblas_LIBRARY)
        include(cmake-modules/FindLAPACKE.cmake)
        set(OpenBLAS_MULTITHREADED 1 )
    endif()
endif()



if(BLAS_openblas_LIBRARY)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "                      ${BLAS_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(blas   INTERFACE IMPORTED)
    add_library(lapack INTERFACE IMPORTED)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${BLAS_openblas_LIBRARY})
    set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
    add_definitions(-DOpenBLAS_AVAILABLE)

else()
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")
    set(OpenBLAS_MULTITHREADED 1 )
    if(OpenMP_FOUND)
        set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    else()
        set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    endif()

    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.3.4
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""

            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) TARGET=${OPENBLAS_MARCH} USE_THREAD=${OpenBLAS_MULTITHREADED} USE_OPENMP=${OpenBLAS_USE_OPENMP} NO_AFFINITY=1 GEMM_MULTITHREAD_THRESHOLD=50 NUM_THREADS=64 BINARY64=64 QUIET_MAKE=1 FFLAGS=-Wno-maybe-uninitialized
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )

ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    add_library(blas   INTERFACE)
    add_library(lapack INTERFACE)
    add_dependencies(blas         library_OpenBLAS)
    add_dependencies(lapack       library_OpenBLAS)
    set(BLAS_INCLUDE_DIRS   ${INSTALL_DIR}/include)
    set(LAPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)

    set(BLAS_LIBRARIES           ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACK_LIBRARIES         ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(BLAS_LIBRARIES_STATIC    ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_STATIC  ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_SHARED    ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_SHARED  ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACKE_FROM_OPENBLAS 1)

endif()

set_target_properties(blas PROPERTIES
        INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES};${PTHREAD_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES   "${BLAS_INCLUDE_DIRS}"
#        INTERFACE_LINK_OPTIONS          "${PTHREAD_LIBRARY}"
        )

set_target_properties(lapack PROPERTIES
        INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES};${PTHREAD_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES   "${LAPACK_INCLUDE_DIRS}"
#        INTERFACE_LINK_OPTIONS          "${PTHREAD_LIBRARY}"
        )


set(FC_LDLAGS -fPIC ${PTHREAD_LIBRARY})

if(OpenBLAS_USE_OPENMP AND OpenMP_FOUND)
    set_target_properties(blas   PROPERTIES INTERFACE_LINK_OPTIONS "${OpenMP_CXX_FLAGS}")
    set_target_properties(lapack PROPERTIES INTERFACE_LINK_OPTIONS "${OpenMP_CXX_FLAGS}")
    list(APPEND BLAS_LIBRARIES        ${OpenMP_LIBRARIES})
    list(APPEND BLAS_LIBRARIES_STATIC ${OpenMP_LIBRARIES})
    list(APPEND BLAS_LIBRARIES_SHARED ${OpenMP_LIBRARIES})
endif()

add_definitions(-DOpenBLAS_AVAILABLE)

