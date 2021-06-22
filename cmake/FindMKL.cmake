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



function(register_found_components)
    foreach(cmp ${MKL_FIND_COMPONENTS})
        if(MKL_${cmp}_FOUND)
            set(MKL_${cmp}_FOUND TRUE PARENT_SCOPE)
        endif()
    endforeach()
endfunction()


function(find_mkl_libraries)

    if(UNIX AND NOT APPLE)
        if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
            set(MKL_ARCH_DIR "intel64")
        else()
            set(MKL_ARCH_DIR "ia32")
        endif()
    endif()

    if(NOT BUILD_SHARED_LIBS)
        # Prefer static libraries
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
    if(BUILD_SHARED_LIBS)
        set(LINK_TYPE SHARED)
        set(MKL_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
        set(LINK_TYPE STATIC)
        set(MKL_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()

    # macos
    if(APPLE)
        set(MKL_ARCH_DIR "intel64")
    endif()

    if(FORCE_BUILD_32BITS)
        set(MKL_ARCH_DIR "ia32")
    endif()

    if (WIN32)
        if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
            set(MKL_ARCH_DIR "intel64")
        else()
            set(MKL_ARCH_DIR "ia32")
        endif()
    endif()


    set(MKL_ROOT_SEARCH_PATHS
            $ENV{MKL_DIR}  ${MKL_DIR}
            $ENV{MKLDIR}   ${MKLDIR}
            $ENV{MKLROOT}  ${MKLROOT}
            $ENV{MKL_ROOT} ${MKL_ROOT}
            $ENV{mkl_root} ${mkl_root}
            $ENV{HOME}/intel/mkl
            /opt/intel/mkl
            /opt/intel
            /opt
            /usr/lib/x86_64-linux-gnu
            /usr
            /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
            "Program Files (x86)/Intel/ComposerXE-2011/mkl"
            )

    set(MKL_PATH_SUFFIXES
            mkl
            intel/mkl
            intel
            )

    find_path(MKL_ROOT_DIR
            include/mkl.h
            PATHS ${MKL_ROOT_SEARCH_PATHS}
            PATH_SUFFIXES
            ${MKL_PATH_SUFFIXES}
            )
    if(MKL_ROOT_DIR)
        message(STATUS "Found MKL ROOT: ${MKL_ROOT_DIR}")
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

        set(par_libnames tbb tbbmalloc iomp5)
        set(mkl_libnames rt core sequential intel_thread gnu_thread tbb_thread)
        if(MKL_ARCH_DIR MATCHES "32")
            list(APPEND mkl_libnames
                    intel
                    gf
                    blas
                    lapack)
            else()
            list(APPEND mkl_libnames
                    intel_lp
                    intel_ilp64
                    gf_lp64
                    gf_ilp64
                    blas95_lp64
                    blas95_ilp64
                    lapack95_lp64
                    lapack95_ilp64
                    )
        endif()

        foreach(lib ${mkl_libnames})
            find_library(MKL_${lib}_LIBRARY
                    mkl_${lib}
                    HINTS ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR}
                    )
            if(MKL_${lib}_LIBRARY)
                add_library(mkl::mkl_${lib} UNKNOWN IMPORTED)
                set_target_properties(mkl::mkl_${lib} PROPERTIES IMPORTED_LOCATION "${MKL_${lib}_LIBRARY}")
                set_target_properties(mkl::mkl_${lib} PROPERTIES LINK_WHAT_YOU_USE TRUE)
#                target_include_directories(mkl::mkl_${lib} SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
            endif()
        endforeach()
        if(TARGET mkl::mkl_rt)
            set_target_properties(mkl::mkl_rt PROPERTIES INTERFACE_LINK_DIRECTORIES  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
        endif()

        if("${MKL_ARCH_DIR}" STREQUAL "32")
            set(TBB_SEARCH ia32)
        else()
            set(TBB_SEARCH intel64)
        endif()
        foreach(lib ${par_libnames})
            find_library(MKL_${lib}_LIBRARY
                    ${lib}
                    HINTS
                    ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR}
                    ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.7
                    ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.4
                    ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR} ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR}
                    )
            if(MKL_${lib}_LIBRARY)
                add_library(mkl::${lib} UNKNOWN IMPORTED)
                set_target_properties(mkl::${lib} PROPERTIES IMPORTED_LOCATION "${MKL_${lib}_LIBRARY}")
                set_target_properties(mkl::${lib} PROPERTIES LINK_WHAT_YOU_USE TRUE)
            endif()
        endforeach()

        # Define usable targets
        set (MKL_FORTRAN_VARIANTS intel gf)
        set (MKL_ARCH_VARIANTS)
        set (MKL_BLAS_SUFFIX)
        set (MKL_THREAD_VARIANTS sequential intel_thread gnu_thread tbb_thread )
        if(MKL_ARCH_DIR MATCHES "64")
            list(APPEND MKL_ARCH_VARIANTS _ilp64 _lp64)
            list(APPEND MKL_BLAS_SUFFIX 95)
        endif()

        foreach (FORTRANVAR ${MKL_FORTRAN_VARIANTS})
            foreach (ARCHVAR ${MKL_ARCH_VARIANTS})
                foreach (THREADVAR ${MKL_THREAD_VARIANTS})
                    foreach(SFX ${MKL_BLAS_SUFFIX})
                        if(THREADVAR MATCHES "sequential")
                            set(threadv seq)
                        elseif(THREADVAR MATCHES "intel_thread")
                            set(threadv ithread)
                        elseif(THREADVAR MATCHES "gnu_thread")
                            set(threadv gthread)
                        elseif(THREADVAR MATCHES "tbb_thread")
                            set(threadv tthread)
                        endif()
                        add_library(mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR} INTERFACE IMPORTED)
                        target_link_libraries(mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR} INTERFACE
                                -Wl,--no-as-needed
                                mkl::mkl_blas${SFX}${ARCHVAR}
                                mkl::mkl_lapack${SFX}${ARCHVAR}
                                -Wl,--start-group
                                mkl::mkl_${FORTRANVAR}${ARCHVAR}
                                mkl::mkl_${THREADVAR}
                                mkl::mkl_core
                                -Wl,--end-group
                                -Wl,--as-needed
                                m
                                dl
                                )
                        target_include_directories(mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR} SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
                        set_target_properties(mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR} PROPERTIES INTERFACE_LINK_DIRECTORIES ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
                        if(MKL_ARCH_DIR MATCHES "64")
                            target_compile_options(mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR} INTERFACE -m64)
                        endif()
                        list(APPEND MKL_TARGETS mkl::mkl_${FORTRANVAR}_${threadv}${ARCHVAR})
                    endforeach()
                endforeach()
            endforeach()
        endforeach()
        set(MKL_TARGETS ${MKL_TARGETS} PARENT_SCOPE)
        set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} PARENT_SCOPE)
    endif()
endfunction()

find_mkl_libraries()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_TARGETS)

