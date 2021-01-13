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

    if(BUILD_SHARED_LIBS)
        set(LINK_TYPE SHARED)
        set(MKL_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
        set(LINK_TYPE STATIC)
        set(MKL_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()


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

        find_library(MKL_CORE_LIBRARY
                mkl_core
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_core ${LINK_TYPE} IMPORTED)
        target_include_directories(mkl::mkl_core SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
        set_target_properties(mkl::mkl_core PROPERTIES IMPORTED_LOCATION "${MKL_CORE_LIBRARY}")
        set_target_properties(mkl::mkl_core PROPERTIES LINK_WHAT_YOU_USE TRUE)

        # Threading libraries
        find_library(MKL_SEQUENTIAL_LIBRARY
                mkl_sequential
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_sequential ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_sequential PROPERTIES IMPORTED_LOCATION "${MKL_SEQUENTIAL_LIBRARY}")
        set_target_properties(mkl::mkl_sequential PROPERTIES LINK_WHAT_YOU_USE TRUE)

        find_library(MKL_INTELTHREAD_LIBRARY
                mkl_intel_thread
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_intel_thread ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_intel_thread PROPERTIES IMPORTED_LOCATION "${MKL_INTELTHREAD_LIBRARY}")
        set_target_properties(mkl::mkl_intel_thread PROPERTIES LINK_WHAT_YOU_USE TRUE)

        find_library(MKL_GNUTHREAD_LIBRARY
                mkl_gnu_thread
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_gnu_thread ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_gnu_thread PROPERTIES IMPORTED_LOCATION "${MKL_GNUTHREAD_LIBRARY}")
        set_target_properties(mkl::mkl_gnu_thread PROPERTIES LINK_WHAT_YOU_USE TRUE)

        find_library(MKL_TBBTHREAD_LIBRARY
                mkl_tbb_thread
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_tbb_thread ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_tbb_thread PROPERTIES IMPORTED_LOCATION "${MKL_TBBTHREAD_LIBRARY}")
        set_target_properties(mkl::mkl_tbb_thread PROPERTIES LINK_WHAT_YOU_USE TRUE)



        # Intel Fortran Libraries
        if(MKL_ARCH_DIR MATCHES "32")
            find_library(MKL_INTEL_LP_LIBRARY
                    mkl_intel
                    HINTS ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR}
                    )

            # Dummy... it's actually the same library as above
            find_library(MKL_INTEL_ILP_LIBRARY
                    mkl_intel${MKL_LIB_SUFFIX}
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
                    mkl_intel_ilp64${MKL_LIB_SUFFIX}
                    HINTS ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR}
                    )

        endif()
        add_library(mkl::mkl_intel_lp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_intel_lp PROPERTIES IMPORTED_LOCATION "${MKL_INTEL_LP_LIBRARY}")
        set_target_properties(mkl::mkl_intel_lp PROPERTIES LINK_WHAT_YOU_USE TRUE)

        add_library(mkl::mkl_intel_ilp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_intel_ilp PROPERTIES IMPORTED_LOCATION "${MKL_INTEL_ILP_LIBRARY}")
        set_target_properties(mkl::mkl_intel_ilp PROPERTIES LINK_WHAT_YOU_USE TRUE)




        # GNU Fortran Libraries
        if(MKL_ARCH_DIR MATCHES "32")
            find_library(MKL_GF_LP_LIBRARY
                    mkl_gf
                    HINTS ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR}
                    )
            find_library(MKL_GF_ILP_LIBRARY
                    mkl_gf
                    HINTS ${MKL_ROOT_DIR}
                    PATH_SUFFIXES
                    lib lib/${MKL_ARCH_DIR}
                    )
        else()
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
        endif()
        add_library(mkl::mkl_gf_lp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_gf_lp PROPERTIES IMPORTED_LOCATION "${MKL_GF_LP_LIBRARY}")
        set_target_properties(mkl::mkl_gf_lp PROPERTIES LINK_WHAT_YOU_USE TRUE)


        add_library(mkl::mkl_gf_ilp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_gf_ilp PROPERTIES IMPORTED_LOCATION "${MKL_GF_ILP_LIBRARY}")
        set_target_properties(mkl::mkl_gf_ilp PROPERTIES LINK_WHAT_YOU_USE TRUE)


        # Blas
        if("${MKL_ARCH_DIR}" STREQUAL "32")
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

        else()
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

        endif()

        add_library(mkl::mkl_blas_lp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_blas_lp PROPERTIES IMPORTED_LOCATION "${MKL_BLAS_LP_LIBRARY}")
        set_target_properties(mkl::mkl_blas_lp PROPERTIES LINK_WHAT_YOU_USE TRUE)

        add_library(mkl::mkl_blas_ilp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_blas_ilp PROPERTIES IMPORTED_LOCATION "${MKL_BLAS_ILP_LIBRARY}")
        set_target_properties(mkl::mkl_blas_ilp PROPERTIES LINK_WHAT_YOU_USE TRUE)



        # Lapack
        if(MKL_ARCH_DIR MATCHES "32")
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
        else()
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

        endif()
        add_library(mkl::mkl_lapack_lp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_lapack_lp PROPERTIES IMPORTED_LOCATION "${MKL_LAPACK_LP_LIBRARY}")
        set_target_properties(mkl::mkl_lapack_lp PROPERTIES LINK_WHAT_YOU_USE TRUE)

        add_library(mkl::mkl_lapack_ilp ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_lapack_ilp PROPERTIES IMPORTED_LOCATION "${MKL_LAPACK_ILP_LIBRARY}")
        set_target_properties(mkl::mkl_lapack_ilp PROPERTIES LINK_WHAT_YOU_USE TRUE)



        # Single shared library
        find_library(MKL_RT_LIBRARY
                mkl_rt
                HINTS ${MKL_ROOT_DIR}
                PATH_SUFFIXES
                lib lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::mkl_rt ${LINK_TYPE} IMPORTED)
        set_target_properties(mkl::mkl_rt PROPERTIES IMPORTED_LOCATION "${MKL_RT_LIBRARY}")
        set_target_properties(mkl::mkl_rt PROPERTIES INTERFACE_LINK_DIRECTORIES  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
        set_target_properties(mkl::mkl_rt PROPERTIES LINK_WHAT_YOU_USE TRUE)
        target_include_directories(mkl::mkl_rt SYSTEM INTERFACE ${MKL_INCLUDE_DIR})


        # iomp5
        if(UNIX AND NOT APPLE)
            find_library(MKL_IOMP5_LIBRARY
                    iomp5
                    HINTS
                    ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR}
                    )
            add_library(mkl::iomp5 ${LINK_TYPE} IMPORTED)
            set_target_properties(mkl::iomp5 PROPERTIES IMPORTED_LOCATION "${MKL_IOMP5_LIBRARY}")
        else()
            SET(MKL_IOMP5_LIBRARY "") # no need for mac
            add_library(mkl::iomp5 INTERFACE IMPORTED)
        endif()

        # TBB
        if("${MKL_ARCH_DIR}" STREQUAL "32")
            set(TBB_SEARCH ia32)
        else()
            set(TBB_SEARCH intel64)
        endif()
        find_library(MKL_TBB_LIBRARY
                NAMES
                tbb libtbb libtbb.so libtbb.so.2
                HINTS
                ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.7
                ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.4
                ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::tbb SHARED IMPORTED)
        set_target_properties(mkl::tbb PROPERTIES IMPORTED_LOCATION "${MKL_TBB_LIBRARY}")
        set_target_properties(mkl::tbb PROPERTIES LINK_WHAT_YOU_USE TRUE)
        # TBBMALLOC
        find_library(MKL_TBBMALLOC_LIBRARY
                NAMES
                tbbmalloc libtbbmalloc libtbbmalloc.so libtbbmalloc.so.2
                HINTS
                ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.7
                ${MKL_ROOT_DIR}/../tbb/lib/${TBB_SEARCH}/gcc4.4
                ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR}
                )
        add_library(mkl::tbbmalloc SHARED IMPORTED)
        set_target_properties(mkl::tbbmalloc PROPERTIES IMPORTED_LOCATION "${MKL_TBBMALLOC_LIBRARY}")
        set_target_properties(mkl::tbbmalloc PROPERTIES LINK_WHAT_YOU_USE TRUE)
        target_link_libraries(mkl::tbb INTERFACE mkl::tbbmalloc)

        # Define usable targets
        set (MKL_FORTRAN_VARIANTS intel gf)
        set (MKL_MODE_VARIANTS ilp lp)
        set (MKL_THREAD_VARIANTS sequential intel_thread gnu_thread tbb_thread )

        foreach (FORTRANVAR ${MKL_FORTRAN_VARIANTS})
            foreach (MODEVAR ${MKL_MODE_VARIANTS})
                foreach (THREADVAR ${MKL_THREAD_VARIANTS})
                    if(THREADVAR MATCHES "sequential")
                        set(threadv seq)
                    endif()
                    if(THREADVAR MATCHES "intel_thread")
                        set(threadv ithread)
                    endif()
                    if(THREADVAR MATCHES "gnu_thread")
                        set(threadv gthread)
                    endif()
                    if(THREADVAR MATCHES "tbb_thread")
                        set(threadv tthread)
                    endif()
                    add_library(mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv} INTERFACE IMPORTED)
                    target_link_libraries(mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv} INTERFACE
#                            -Wl,--no-as-needed
                            mkl::mkl_blas_${MODEVAR}
                            mkl::mkl_lapack_${MODEVAR}
                            -Wl,--start-group
                            mkl::mkl_${FORTRANVAR}_${MODEVAR}
                            mkl::mkl_${THREADVAR}
                            mkl::mkl_core
                            -Wl,--end-group
#                            -Wl,--as-needed
                            m
                            dl
                            )
                    target_include_directories(mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv} SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
                    set_target_properties(mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv} PROPERTIES INTERFACE_LINK_DIRECTORIES  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
                    if(MKL_ARCH_DIR MATCHES "intel64")
                        target_compile_options(mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv} INTERFACE -m64)
                    endif()
                    list(APPEND MKL_TARGETS mkl::mkl_${FORTRANVAR}_${MODEVAR}_${threadv})
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

