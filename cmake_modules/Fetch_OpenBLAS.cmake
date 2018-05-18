
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

# It seems that libopenblas from "apt" includes lapack in newer versions of Ubuntu. Trusty LTS 14.04 does not
# libopenblas bundled with lapack. Therefore we test if both lapack and blas are present to distinghuish these cases.
# Otherwise, arpack-ng will complain later about undefined references.

enable_language(Fortran)
include(cmake_modules/FindGFortran.cmake)

message(STATUS "SEARCHING FOR OpenBLAS IN SYSTEM...")
set(BLA_VENDOR Open)
set(BLAS_VERBOSE OFF)
find_package(BLAS)
set(BLA_VENDOR OpenBLAS)
find_package(LAPACK)
if(BLAS_FOUND AND LAPACK_FOUND)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "LAPACK FOUND IN SYSTEM: ${LAPACK_openblas_LIBRARY}")
       #For convenience, define these variables
    add_library(blas INTERFACE)
    add_library(lapack INTERFACE)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${LAPACK_openblas_LIBRARY})
endif()
#exit (1)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.2.20
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) USE_THREAD=0 USE_OPENMP=0 NO_LAPACKE=1 NO_CBLAS=1 BINARY64=1 QUIET_MAKE=1
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            )

    ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    add_library(blas INTERFACE)
    add_library(lapack INTERFACE)
    add_dependencies(blas         library_OpenBLAS)
    add_dependencies(lapack       library_OpenBLAS)
    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(BLAS_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})

endif()

set_target_properties(blas PROPERTIES
        INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            ""
        INTERFACE_COMPILE_OPTIONS       ""
        )

set_target_properties(lapack PROPERTIES
        INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            ""
        INTERFACE_COMPILE_OPTIONS       ""
        )



target_link_libraries(${PROJECT_NAME} PUBLIC blas)
target_link_libraries(${PROJECT_NAME} PUBLIC lapack)
target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})