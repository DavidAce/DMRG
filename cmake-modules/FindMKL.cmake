# - Try to find the Intel Math Kernel Library
#   Forked from:  https://github.com/Eyescale/CMake/blob/master/FindMKL.cmake
#   which is forked from: https://github.com/openmeeg/openmeeg/blob/master/macros/FindMKL.cmake

# Once done this will define
#
# MKL_FOUND - system has MKL
# MKL_ROOT_DIR - path to the MKL base directory
# MKL_INCLUDE_DIR - the MKL include directory
# MKL_LIBRARIES - MKL libraries
#
# There are few sets of libraries:
# Fortran modes:
# GF - GNU Fortran
# INTEL - INTEL Fortran
#
# Array indexes modes:
# LP - 32 bit indexes of arrays
# ILP - 64 bit indexes of arrays
#
# Threading:
# SEQUENTIAL - no threading
# INTEL - Intel threading library
# GNU - GNU threading library
#
# MPI support
# NOMPI - no MPI support
# INTEL - Intel MPI library
# OPEN - Open MPI library
# SGI - SGI MPT Library

# linux
if(UNIX AND NOT APPLE)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        set(MKL_ARCH_DIR "intel64")
    else()
        set(MKL_ARCH_DIR "32")
    endif()
endif()

# macos
if(APPLE)
    set(MKL_ARCH_DIR "intel64")
endif()

IF(FORCE_BUILD_32BITS)
    set(MKL_ARCH_DIR "32")
ENDIF()

if (WIN32)
    if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
        set(MKL_ARCH_DIR "intel64")
    else()
        set(MKL_ARCH_DIR "ia32")
    endif()
endif()

set (MKL_FORTRAN_VARIANTS INTEL GF)
set (MKL_MODE_VARIANTS ILP LP)
set (MKL_THREAD_VARIANTS SEQUENTIAL GNUTHREAD INTELTHREAD)
set (MKL_MPI_VARIANTS NOMPI INTELMPI OPENMPI SGIMPT)

set(MKL_ROOT_SEARCH_PATHS
        $ENV{MKL_DIR}  ${MKL_DIR}
        $ENV{MKLDIR}   ${MKLDIR}
        $ENV{MKLROOT}  ${MKLROOT}
        $ENV{MKL_ROOT} ${MKL_ROOT}
        $ENV{mkl_root} ${mkl_root}
        $ENV{HOME}/intel/mkl
        /opt/intel/mkl
        $ENV{BLAS_DIR}
        /usr/lib/x86_64-linux-gnu
        /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
        "Program Files (x86)/Intel/ComposerXE-2011/mkl"
        )
set(MKL_PATH_SUFFIXES
        mkl
        intel/mkl

        )



find_path(MKL_ROOT_DIR
        include/mkl.h
        HINTS ${DIRECTORY_HINTS}
        PATHS ${MKL_ROOT_SEARCH_PATHS}
        )

find_path(MKL_INCLUDE_DIR
        mkl.h
        HINTS ${MKL_ROOT_DIR}/include
        )

find_path(MKL_FFTW_INCLUDE_DIR
        fftw3.h
        HINTS ${MKL_ROOT_DIR}/include
        PATH_SUFFIXES fftw
        NO_DEFAULT_PATH
        )

find_library(MKL_CORE_LIBRARY
        mkl_core
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )

# Threading libraries
find_library(MKL_SEQUENTIAL_LIBRARY
        mkl_sequential
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )

find_library(MKL_INTELTHREAD_LIBRARY
        mkl_intel_thread
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )

find_library(MKL_GNUTHREAD_LIBRARY
        mkl_gnu_thread
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )

# Intel Fortran Libraries
IF("${MKL_ARCH_DIR}" STREQUAL "32")
    find_library(MKL_INTEL_LP_LIBRARY
            mkl_intel
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )

    find_library(MKL_INTEL_ILP_LIBRARY
            mkl_intel
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
else()
    find_library(MKL_INTEL_LP_LIBRARY
            mkl_intel_lp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )

    find_library(MKL_INTEL_ILP_LIBRARY
            mkl_intel_ilp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
ENDIF()

# GNU Fortran Libraries
find_library(MKL_GF_LP_LIBRARY
        mkl_gf_lp64
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )

find_library(MKL_GF_ILP_LIBRARY
        mkl_gf_ilp64
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )



# Blas
IF("${MKL_ARCH_DIR}" STREQUAL "32")
    find_library(MKL_BLAS_LP_LIBRARY
            mkl_blas
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
    find_library(MKL_BLAS_ILP_LIBRARY
            mkl_blas
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
ELSE()
    find_library(MKL_BLAS_LP_LIBRARY
            mkl_blas95_lp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
    find_library(MKL_BLAS_ILP_LIBRARY
            mkl_blas95_ilp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )

ENDIF()



# Lapack
IF("${MKL_ARCH_DIR}" STREQUAL "32")
    find_library(MKL_LAPACK_LP_LIBRARY
            mkl_lapack
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
    find_library(MKL_LAPACK_ILP_LIBRARY
            mkl_lapack
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
ELSE()
    find_library(MKL_LAPACK_LP_LIBRARY
            mkl_lapack95_lp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )
    find_library(MKL_LAPACK_ILP_LIBRARY
            mkl_lapack95_ilp64
            HINTS ${MKL_ROOT_DIR}
            PATH_SUFFIXES
            lib lib/${MKL_ARCH_DIR}
            )

ENDIF()


# Single shared library
find_library(MKL_RT_LIBRARY
        mkl_rt
        HINTS ${MKL_ROOT_DIR}
        PATH_SUFFIXES
        lib lib/${MKL_ARCH_DIR}
        )


# iomp5
IF("${MKL_ARCH_DIR}" STREQUAL "32")
    IF(UNIX AND NOT APPLE)
        find_library(MKL_IOMP5_LIBRARY
                iomp5
                PATHS
                ${MKL_ROOT_DIR}/../lib/ia32
                )
    ELSE()
        SET(MKL_IOMP5_LIBRARY "") # no need for mac
    ENDIF()
else()
    IF(UNIX AND NOT APPLE)
        find_library(MKL_IOMP5_LIBRARY
                iomp5
                PATHS
                ${MKL_ROOT_DIR}/../lib/intel64
                )
    ELSE()
        SET(MKL_IOMP5_LIBRARY "") # no need for mac
    ENDIF()
ENDIF()

foreach (FORTRANVAR ${MKL_FORTRAN_VARIANTS})
    foreach (MODEVAR ${MKL_MODE_VARIANTS})
        foreach (THREADVAR ${MKL_THREAD_VARIANTS})
            if (MKL_CORE_LIBRARY AND MKL_${FORTRANVAR}_${MODEVAR}_LIBRARY AND MKL_${THREADVAR}_LIBRARY)
                set(MKL_${FORTRANVAR}_${MODEVAR}_${THREADVAR}_LIBRARIES
                        ${MKL_${FORTRANVAR}_${MODEVAR}_LIBRARY}
#                        ${MKL_${MODEVAR}_LIBRARY}
                        ${MKL_${THREADVAR}_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_LAPACK_LIBRARY}
                        ${MKL_IOMP5_LIBRARY})
#                message("${FORTRANVAR} ${MODEVAR} ${THREADVAR} ${MKL_${FORTRANVAR}_${MODEVAR}_${THREADVAR}_LIBRARIES}") # for debug
            endif()
        endforeach()
    endforeach()
endforeach()

set(MKL_LIBRARIES ${MKL_GF_LP_GNUTHREAD_LIBRARIES})

if($ENV{LD_LIBRARY_PATH} MATCHES "mkl")
    message(STATUS "Intel MKL library path found in LD_LIBRARY_PATH")
else()
    message(STATUS "Intel MKL library path not found your in LD_LIBRARY_PATH. Make sure the path is known by ld when linking.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES
        MKL_CORE_LIBRARY MKL_INTEL_LP_LIBRARY MKL_INTEL_ILP_LIBRARY MKL_GF_LP_LIBRARY MKL_GF_ILP_LIBRARY
        MKL_SEQUENTIAL_LIBRARY MKL_INTELTHREAD_LIBRARY MKL_GNUTHREAD_LIBRARY
        )