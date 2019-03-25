
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

# It seems that libopenblas from "apt" includes lapack in newer versions of Ubuntu. Trusty LTS 14.04 does not
# libopenblas bundled with lapack. Therefore we test if both lapack and blas are present to distinghuish these cases.
# Otherwise, arpack-ng will complain later about undefined references.


set(BLA_VENDOR OpenBLAS)
set(BLAS_VERBOSE OFF)
set(BLA_STATIC ${STATIC_BUILD})

if(NOT BLAS_FOUND)
    message(STATUS "Searching for OpenBLAS in system.")
    find_package(OpenBLAS 0.3
            PATHS
            ${INSTALL_DIRECTORY}/OpenBLAS
            $ENV{BLAS_DIR}/lib
            $ENV{HOME}/.conda/lib
            $ENV{HOME}/.conda
            $ENV{HOME}/anaconda3/lib
            $ENV{HOME}/anaconda3
            /usr/lib/x86_64-linux-gnu/
            )
    if(OpenBLAS_FOUND)
        get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
        get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
        set(BLAS_openblas_LIBRARY ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
        set(BLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
    endif()
endif()

if(NOT BLAS_openblas_LIBRARY)
    message(STATUS "Searching for OpenBLAS in standard paths")
    find_library(BLAS_openblas_LIBRARY
            NAMES libopenblas${CUSTOM_SUFFIX}
            PATHS
            ${INSTALL_DIRECTORY}/OpenBLAS
            $ENV{BLAS_DIR}/lib
            $ENV{HOME}/.conda/lib
            $ENV{HOME}/anaconda3/lib
            /usr/lib/x86_64-linux-gnu
            NO_DEFAULT_PATH
            )
    find_path(BLAS_INCLUDE_DIRS
            NAMES cblas.h
            PATHS
            ${INSTALL_DIRECTORY}/OpenBLAS
            $ENV{BLAS_DIR}/include
            $ENV{HOME}/.conda/include
            $ENV{HOME}/anaconda3/include
            /usr/include
            /usr/include/x86_64-linux-gnu
            NO_DEFAULT_PATH
            )
    if (BLAS_openblas_LIBRARY)
        set(OpenBLAS_MULTITHREADED 1 )
        set(OpenBLAS_MULTITHREADED 1 )
    endif()
endif()




# To print all variables, use the code below:
#
# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     if("${_variableName}" MATCHES "BLAS" OR "${_variableName}" MATCHES "blas" OR "${_variableName}" MATCHES "Blas")
#         message(STATUS "${_variableName}=${${_variableName}}")
#     endif()
# endforeach()



if(BLAS_openblas_LIBRARY)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "                      ${BLAS_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(blas   INTERFACE)
    add_library(lapack INTERFACE)
    set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
    set(LAPACK_LIBRARIES   ${BLAS_openblas_LIBRARY})
    set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})

    include(CheckIncludeFileCXX)
    include(CheckSymbolExists)
    check_include_file_cxx(cblas.h    has_cblas_h)
    if(has_cblas_h)
        set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${PTHREAD_LIBRARY} gfortran)
        set(CMAKE_REQUIRED_INCLUDES  ${BLAS_INCLUDE_DIRS})
        check_symbol_exists(openblas_set_num_threads cblas.h has_openblas_symbols)
    endif()
    if(has_cblas_h AND has_openblas_symbols)
        add_definitions(-DOpenBLAS_AVAILABLE)
    endif()
else()
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")
    message(STATUS "OpenBLAS TARGET: ${OPENBLAS_MARCH}")
    set(OpenBLAS_MULTITHREADED 1 )
    if(OpenMP_FOUND)
        set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    else()
        set(OpenBLAS_USE_OPENMP 0) # Openmp doesnt work on clang it seems
    endif()

    include(ExternalProject)
    ExternalProject_Add(external_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.3.5
            PREFIX      ${BUILD_DIRECTORY}/OpenBLAS
            INSTALL_DIR ${INSTALL_DIRECTORY}/OpenBLAS
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""

            BUILD_IN_SOURCE 1
            BUILD_COMMAND export LD_LIBRARY_PATH=${GFORTRAN_PATH} &&
                          export LDFLAGS=-L${GFORTRAN_LIB_SHARED} &&
                          $(MAKE) TARGET=${OPENBLAS_MARCH}
                          USE_THREAD=${OpenBLAS_MULTITHREADED}
                          USE_OPENMP=${OpenBLAS_USE_OPENMP}
                          NO_AFFINITY=1
                          NO_WARMUP=1
                          QUIET_MAKE=0
                          NUM_THREADS=128
#                          FFLAGS=-frecursive
                          DYNAMIC_ARCH=1
                          BINARY64=64
                          GEMM_MULTITHREAD_THRESHOLD=4
                          #LDFLAGS=-L${GFORTRAN_LIB_SHARED}
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )

ExternalProject_Get_Property(external_OpenBLAS INSTALL_DIR)
    add_library(blas   INTERFACE)
    add_library(lapack INTERFACE)
    add_library(openblas::lapacke ALIAS blas)
    add_dependencies(blas         external_OpenBLAS)
    add_dependencies(lapack       external_OpenBLAS)
    set(BLAS_INCLUDE_DIRS   ${INSTALL_DIR}/include)
    set(LAPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)

    set(BLAS_LIBRARIES           ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(LAPACK_LIBRARIES         ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(BLAS_LIBRARIES_STATIC    ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_STATIC  ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BLAS_LIBRARIES_SHARED    ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES_SHARED  ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_definitions(-DOpenBLAS_AVAILABLE)

endif()



target_link_libraries(blas INTERFACE ${BLAS_LIBRARIES} ${PTHREAD_LIBRARY} gfortran)
target_include_directories(blas INTERFACE ${BLAS_INCLUDE_DIRS})

target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES} ${PTHREAD_LIBRARY} gfortran)
target_include_directories(lapack INTERFACE ${LAPACK_INCLUDE_DIRS})

set(FC_LDLAGS -fPIC ${PTHREAD_LIBRARY})

if(OpenBLAS_USE_OPENMP AND OpenMP_FOUND)
    target_link_libraries(blas INTERFACE ${OpenMP_LIBRARIES})
    target_link_libraries(lapack INTERFACE ${OpenMP_LIBRARIES})
    target_link_options(blas INTERFACE ${OpenMP_CXX_FLAGS})
    target_link_options(lapack INTERFACE ${OpenMP_CXX_FLAGS})
endif()




if (TARGET blas AND TARGET lapack)
    include(cmake-modules/FindLAPACKE.cmake)
    if(TARGET lapacke)
        target_link_libraries(blas   INTERFACE lapacke)
        target_link_libraries(lapack INTERFACE lapacke)
    endif()
endif()
