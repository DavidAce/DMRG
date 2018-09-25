
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


if (EXISTS "$ENV{BLAS_DIR}")
    # Try finding armadillo as module library
    message(STATUS "Attempting to find BLAS and LAPACK from environment modules.")
    find_library(BLAS_openblas_LIBRARY
            NAMES libopenblas.a
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
            NAMES libopenblas.a
            PATHS "/usr/lib/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
    #    find_library(LAPACK_openblas_LIBRARY
    #            NAMES liblapacke.a
    #            PATHS "/usr/lib/x86_64-linux-gnu"
    #            NO_DEFAULT_PATH
    #            )
    find_path(BLAS_INCLUDE_DIRS
            NAMES cblas.h
            PATHS "/usr/include" "/usr/include/x86_64-linux-gnu"
            NO_DEFAULT_PATH
            )
    #    find_path(LAPACK_INCLUDE_DIRS
    #            NAMES lapacke.h
    #            PATHS "/usr/include" "/usr/include/x86_64-linux-gnu"
    #            NO_DEFAULT_PATH
    #            )
    if (BLAS_openblas_LIBRARY)
        set(OpenBLAS_MULTITHREADED 1 )
    endif()
endif()



if(BLAS_openblas_LIBRARY)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "                      ${BLAS_INCLUDE_DIRS}")


    #    message(STATUS "LAPACK FOUND IN SYSTEM: ${LAPACK_openblas_LIBRARY}")
    #    message(STATUS "                        ${LAPACK_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(blas STATIC IMPORTED)
    add_library(lapack STATIC IMPORTED)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${BLAS_openblas_LIBRARY})
    set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
    add_definitions(-DOpenBLAS_AVAILABLE)

else()
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")
    set(OpenBLAS_MULTITHREADED 1 )
    set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.3.3
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) USE_THREAD=${OpenBLAS_MULTITHREADED} USE_OPENMP=${OpenBLAS_USE_OPENMP} OPENBLAS_NUM_THREADS=8 NUM_THREADS=8 BINARY64=1 QUIET_MAKE=1 TARGET=${OPENBLAS_MARCH}
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )
    ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas         library_OpenBLAS)
    add_dependencies(lapack       library_OpenBLAS)
    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(LAPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)

    set(BLAS_LIBRARIES           ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES         ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_STATIC    ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_STATIC  ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_SHARED    ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_SHARED  ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACKE_FROM_OPENBLAS 1)

endif()
set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION               "${BLAS_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES};gfortran"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            "-lpthread"
        )

set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION               "${LAPACK_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES};gfortran"
        INTERFACE_INCLUDE_DIRECTORY     "${LAPACK_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            "-lpthread"
        )


target_link_libraries(${PROJECT_NAME} PRIVATE blas lapack )
target_include_directories(${PROJECT_NAME} PRIVATE ${BLAS_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACK_INCLUDE_DIRS})
if(OpenBLAS_MULTITHREADED)
    target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)
#    target_link_libraries(${PROJECT_NAME} PRIVATE -lomp -lpthread )
#    target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp=libomp)
    #    list(APPEND BLAS_LIBRARIES        "-lpthread")
    #    list(APPEND BLAS_LIBRARIES_STATIC "-lpthread")
    #    list(APPEND BLAS_LIBRARIES_SHARED "-lpthread")
#    set(EXTRA_LDLAGS "-lpthread")
    # Fortran linker arguments
    set(FC_LDLAGS "-lpthread")
endif()
add_definitions(-DOpenBLAS_AVAILABLE)

