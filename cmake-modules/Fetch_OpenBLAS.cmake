
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

# It seems that libopenblas from "apt" includes lapack in newer versions of Ubuntu. Trusty LTS 14.04 does not
# libopenblas bundled with lapack. Therefore we test if both lapack and blas are present to distinghuish these cases.
# Otherwise, arpack-ng will complain later about undefined references.


set(BLA_VENDOR OpenBLAS)
set(BLAS_VERBOSE OFF)
if(NOT ${BUILD_SHARED_LIBS})
    set(BLA_STATIC ON)
else()
    set(BLA_STATIC OFF)
endif()

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

    if(OpenBLAS_FOUND AND OpenBLAS_LIBRARIES AND NOT OpenBLAS_LIBRARY)
        get_filename_component(OpenBLAS_LIBRARIES_WE ${OpenBLAS_LIBRARIES} NAME_WE)
        get_filename_component(OpenBLAS_ROOT ${OpenBLAS_LIBRARIES} DIRECTORY)
        set(OpenBLAS_LIBRARY ${OpenBLAS_ROOT}/${OpenBLAS_LIBRARIES_WE}${CUSTOM_SUFFIX})
    endif()

    if(OpenBLAS_FOUND AND OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIRS)
        set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
        set(OpenBLAS_FOUND TRUE)
        set(BLAS_FOUND TRUE)
    endif()
endif()




if(NOT OpenBLAS_FOUND)
    message(STATUS "Searching for OpenBLAS in standard paths")
    find_library(OpenBLAS_LIBRARY
            NAMES libopenblas${CUSTOM_SUFFIX}
            PATHS
            ${INSTALL_DIRECTORY}/OpenBLAS
            $ENV{BLAS_DIR}/lib
            $ENV{HOME}/.conda/lib
            $ENV{HOME}/anaconda3/lib
            /usr/lib/x86_64-linux-gnu
            NO_DEFAULT_PATH
            )
    find_path(OpenBLAS_INCLUDE_DIRS
            NAMES openblas_config.h
            PATHS
            ${INSTALL_DIRECTORY}/OpenBLAS
            $ENV{BLAS_DIR}/include
            $ENV{HOME}/.conda/include
            $ENV{HOME}/anaconda3/include
            /usr/include
            /usr/include/x86_64-linux-gnu
            NO_DEFAULT_PATH
            )
    if (OpenBLAS_LIBRARY AND OpenBLAS_INCLUDE_DIRS)
        set(OpenBLAS_MULTITHREADED 1 )
        set(BLAS_FOUND TRUE)
        set(OpenBLAS_FOUND TRUE)
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



if(OpenBLAS_FOUND)
    message(STATUS "OpenBLAS FOUND IN SYSTEM: ${OpenBLAS_LIBRARY}")
    message(STATUS "                          ${OpenBLAS_INCLUDE_DIRS}")

    #For convenience, define these variables
    add_library(blas   INTERFACE)
    add_library(lapack ALIAS blas)

    include(CheckIncludeFileCXX)
    include(CheckSymbolExists)
    check_include_file_cxx(cblas.h    has_cblas_h)
    if(has_cblas_h)
        set(CMAKE_REQUIRED_LIBRARIES ${OpenBLAS_LIBRARY} ${PTHREAD_LIBRARY} gfortran)
        set(CMAKE_REQUIRED_INCLUDES  ${OpenBLAS_INCLUDE_DIRS})
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

#            CMAKE_ARGS
#            -DTARGET=${OPENBLAS_MARCH}
#            -DDYNAMIC_ARCH:BOOL=OFF
#            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#            -DBINARY=64
#            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
#            -DCMAKE_C_FLAGS=-w
#            -DCMAKE_CXX_FLAGS=-w
#            -DCMAKE_FORTRAN_FLAGS=-w



#            BUILD_COMMAND export LD_LIBRARY_PATH=${GFORTRAN_PATH}:$LD_LIBRARY_PATH &&
#            export LDFLAGS="-L${GFORTRAN_PATH} $LDFLAGS" &&
            BUILD_IN_SOURCE 1
            BUILD_COMMAND
                        export LD_LIBRARY_PATH=${GFORTRAN_PATH} &&
#                        export LDFLAGS=-L${GFORTRAN_PATH} -l${GFORTRAN_LIB} -l${QUADMATH_LIB} &&
                        $(MAKE) TARGET=SANDYBRIDGE
                        DYNAMIC_ARCH=1
                        USE_THREAD=${OpenBLAS_MULTITHREADED}
                        USE_OPENMP=${OpenBLAS_USE_OPENMP}
                        NO_AFFINITY=1
                        NO_WARMUP=1
                        QUIET_MAKE=0
                        NUM_THREADS=128
                        DEBUG=0
                        FFLAGS=-frecursive
                        BINARY64=64
                        GEMM_MULTITHREAD_THRESHOLD=16
                        LDFLAGS=-L${GFORTRAN_PATH} -l${GFORTRAN_LIB} -l${QUADMATH_LIB}
                        FFLAGS=-O3 -Wno-maybe-uninitialized -Wno-conversion -Wno-unused-but-set-variable -Wno-unused-variable
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            DEPENDS gfortran
            )

ExternalProject_Get_Property(external_OpenBLAS INSTALL_DIR)
    add_library(blas   INTERFACE)
    add_library(lapack ALIAS blas)
    add_library(OpenBLAS::lapacke ALIAS blas)
    add_dependencies(blas         external_OpenBLAS)

    set(OpenBLAS_LIBRARY           ${INSTALL_DIR}/lib/libopenblas${CUSTOM_SUFFIX})
    set(OpenBLAS_LIBRARY_STATIC    ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(OpenBLAS_LIBRARY_SHARED    ${INSTALL_DIR}/lib/libopenblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(OpenBLAS_INCLUDE_DIRS      ${INSTALL_DIR}/include)
    add_definitions(-DOpenBLAS_AVAILABLE)
endif()



target_link_libraries(blas INTERFACE ${OpenBLAS_LIBRARY} pthread gfortran)
target_include_directories(blas INTERFACE ${OpenBLAS_INCLUDE_DIRS})
set(FC_LDLAGS ${PTHREAD_LIBRARY} ${GFORTRAN_LIB})

if(OpenBLAS_USE_OPENMP AND OpenMP_FOUND)
    target_link_libraries(blas INTERFACE ${OpenMP_LIBRARIES})
    target_link_options(blas INTERFACE ${OpenMP_CXX_FLAGS})
endif()


if (TARGET blas)
    set(BLAS_LIBRARIES   ${OpenBLAS_LIBRARY} ${FC_LDLAGS})
    set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARY} ${FC_LDLAGS})
    include(cmake-modules/FindLAPACKE.cmake)
    if(TARGET lapacke)
        target_link_libraries(blas   INTERFACE lapacke)
    endif()
endif()
