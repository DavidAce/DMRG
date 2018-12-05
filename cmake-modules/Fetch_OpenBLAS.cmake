
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
    # Try finding armadillo as module library
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

    #For convenience, define these variables
    add_library(blas   UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${BLAS_openblas_LIBRARY})
    set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
    add_definitions(-DOpenBLAS_AVAILABLE)

else()
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")
    set(OpenBLAS_MULTITHREADED 1 )
    set(OpenBLAS_USE_OPENMP 1) # Openmp doesnt work on clang it seems
    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.3.4
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
#            CMAKE_ARGS
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DBUILD_STATIC_LIBS:BOOL=ON
#            #-DBUILD_SHARED_LIBS:BOOL=ON
#            -DUSE_THREAD:BOOL=ON
#            -DUSE_OPENMP:BOOL=ON
#            -DNO_AFFINITY:BOOL=ON
#            -DBINARY64=64
#            -DNUM_THREADS=64
#            -DCMAKE_ANSI_CFLAGS:STRING=-fPIC

            #-DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            #-DCMAKE_C_FLAGS=-w

            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) TARGET=${OPENBLAS_MARCH} USE_THREAD=${OpenBLAS_MULTITHREADED} USE_OPENMP=${OpenBLAS_USE_OPENMP} NO_AFFINITY=1 GEMM_MULTITHREAD_THRESHOLD=50 NUM_THREADS=64 BINARY64=64 QUIET_MAKE=1 CFLAGS=-O3 FFLAGS=-O3
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )

ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    add_library(blas   UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas         library_OpenBLAS)
    add_dependencies(lapack       library_OpenBLAS)
    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(LAPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)

    set(BLAS_LIBRARIES           ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACK_LIBRARIES         ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(BLAS_LIBRARIES_STATIC    ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACK_LIBRARIES_STATIC  ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(BLAS_LIBRARIES_SHARED    ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACK_LIBRARIES_SHARED  ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACKE_FROM_OPENBLAS 1)

endif()
set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION               "${BLAS_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES};gfortran;-lpthread"
        INTERFACE_INCLUDE_DIRECTORY     "${BLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            "-lpthread"
        )

set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION               "${LAPACK_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES};gfortran;-lpthread"
        INTERFACE_INCLUDE_DIRECTORY     "${LAPACK_INCLUDE_DIRS}"
        INTERFACE_LINK_FLAGS            "-lpthread"
        )


target_link_libraries(${PROJECT_NAME} PRIVATE blas lapack )
target_include_directories(${PROJECT_NAME} PRIVATE ${BLAS_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${LAPACK_INCLUDE_DIRS})
if(OpenBLAS_MULTITHREADED)
    target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)
    set_target_properties(blas   PROPERTIES INTERFACE_LINK_LIBRARIES "-lpthread")
    set_target_properties(lapack PROPERTIES INTERFACE_LINK_LIBRARIES "-lpthread")

        #    target_link_libraries(${PROJECT_NAME} PRIVATE -lomp -lpthread )
#    target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp=libomp)
    #    list(APPEND BLAS_LIBRARIES        "-lpthread")
    #    list(APPEND BLAS_LIBRARIES_STATIC "-lpthread")
    #    list(APPEND BLAS_LIBRARIES_SHARED "-lpthread")
#    set(EXTRA_LDLAGS "-lpthread")
    # Fortran linker arguments
    set(FC_LDLAGS "-lpthread")
endif()
if(OpenBLAS_USE_OPENMP AND OpenMP_FOUND)
    set_target_properties(blas   PROPERTIES INTERFACE_LINK_LIBRARIES "${OpenMP_CXX_FLAGS}")
    set_target_properties(lapack PROPERTIES INTERFACE_LINK_LIBRARIES "${OpenMP_CXX_FLAGS}")
endif()


add_definitions(-DOpenBLAS_AVAILABLE)

