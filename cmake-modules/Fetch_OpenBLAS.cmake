
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
set(BLA_STATIC ON)
find_package(BLAS)
find_package(LAPACK)
#set(BLA_VENDOR OpenBLAS)
#find_package(LAPACK)

if (NOT BLAS_openblas_LIBRARY OR NOT LAPACK_openblas_LIBRARY)
    # Try finding armadillo as module library
    find_library(BLAS_openblas_LIBRARY
            NAMES libopenblas.a
            PATHS "$ENV{BLAS_DIR}/lib"
            NO_DEFAULT_PATH
            )
    find_path(BLAS_INCLUDE_DIRS
            NAMES cblas.h
            PATHS "$ENV{BLAS_DIR}/include"
            NO_DEFAULT_PATH
            )

    find_library(LAPACK_openblas_LIBRARY
            NAMES libopenblas.a
            PATHS "$ENV{BLAS_DIR}/lib"
            NO_DEFAULT_PATH
            )
    find_path(LAPACK_INCLUDE_DIRS
            NAMES cblas.h
            PATHS "$ENV{BLAS_DIR}/include"
            NO_DEFAULT_PATH
            )
    if (BLAS_openblas_LIBRARY AND LAPACK_openblas_LIBRARY)
        set(BLAS_FOUND 1)
        set(LAPACK_FOUND 1)
    endif()
endif()





if(BLAS_openblas_LIBRARY AND LAPACK_openblas_LIBRARY)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "                      ${BLAS_INCLUDE_DIRS}")
    message(STATUS "LAPACK FOUND IN SYSTEM: ${LAPACK_openblas_LIBRARY}")
    message(STATUS "                        ${LAPACK_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(blas STATIC IMPORTED)
    add_library(lapack STATIC IMPORTED)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${LAPACK_openblas_LIBRARY})
    add_definitions(-DOpenBLAS_AVAILABLE)

endif()
#exit (1)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")
    set(OpenBLAS_MULTITHREADED 0 )
    set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.3.0
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) USE_THREAD=${OpenBLAS_MULTITHREADED} USE_OPENMP=${OpenBLAS_USE_OPENMP} OPENBLAS_NUM_THREADS=1 NUM_THREADS=1 BINARY64=1 QUIET_MAKE=1 TARGET=NEHALEM
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )
    #NO_LAPACKE=${OpenBLAS_USE_OTHER} NO_CBLAS=${OpenBLAS_USE_OTHER} BINARY64=1 QUIET_MAKE=0
    ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas         library_OpenBLAS)
    add_dependencies(lapack       library_OpenBLAS)
    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(BLAS_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_STATIC ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_STATIC ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_SHARED ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_SHARED ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()
set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION               "${BLAS_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES};gfortran"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            ""
        INTERFACE_COMPILE_OPTIONS       ""
        )

set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION               "${LAPACK_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES};gfortran"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            ""
        INTERFACE_COMPILE_OPTIONS       ""
        )



target_link_libraries(${PROJECT_NAME} PRIVATE blas lapack)
target_include_directories(${PROJECT_NAME} PRIVATE ${BLAS_INCLUDE_DIRS})
add_definitions(-DOpenBLAS_AVAILABLE)

if(OpenBLAS_MULTITHREADED)
    target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)
endif()