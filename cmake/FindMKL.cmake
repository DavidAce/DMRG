# find_package module for the Intel Math Kernel Library (MKL)
#
# COMPONENTS
# architecture:

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
    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES .so.1 .so.2) # For tbb


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
    set(MKL_ARCH_DIR ${MKL_ARCH_DIR} PARENT_SCOPE)

    set(MKL_ROOT_SEARCH_PATHS
            ${MKL_ROOT_DIR}
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

        set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} PARENT_SCOPE)
        if(NOT MKL_INCLUDE_DIR)
            message(WARNING "Found MKL_ROOT_DIR but not MKL_INCLUDE_DIR:"
                    "MKL_ROOT_DIR   : ${MKL_ROOT_DIR})"
                    "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR})")
                    return()
        endif()
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
                    intel_lp64
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
#                message(STATUS "Added library mkl::mkl_${lib}")
            else()
#                message(STATUS "Failed library mkl::mkl_${lib}")

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
                    ${MKL_ROOT_DIR}
                    ${MKL_ROOT_DIR}/../tbb
                    PATH_SUFFIXES
                    lib  ../lib/${MKL_ARCH_DIR}
                    lib/${MKL_ARCH_DIR}
                    lib/${MKL_ARCH_DIR}/gcc4.7
                    lib/${MKL_ARCH_DIR}/gcc4.4
                    )
            if(MKL_${lib}_LIBRARY)
                add_library(mkl::${lib} UNKNOWN IMPORTED)
                set_target_properties(mkl::${lib} PROPERTIES IMPORTED_LOCATION "${MKL_${lib}_LIBRARY}")
                set_target_properties(mkl::${lib} PROPERTIES LINK_WHAT_YOU_USE TRUE)
#                message(STATUS "Added library mkl::${lib}")
            else()
#                message(FATAL_ERROR "Failed library mkl::${lib}")
            endif()
        endforeach()
    endif()
endfunction()

function(setup_mkl_targets)

    # Setup all the library variants
    set (MKL_FORTRAN_VARIANTS intel gf)
    set (MKL_THREAD_VARIANTS sequential intel_thread gnu_thread tbb_thread )
    set (MKL_ARCH_VARIANTS)
    set (MKL_MATH_VARIANTS blas lapack)
    if(MKL_ARCH_DIR MATCHES "64")
        list(APPEND MKL_ARCH_VARIANTS ilp64 lp64)
        list(APPEND MKL_95_SUFFIX 95_)
    else()
        list(APPEND MKL_ARCH_VARIANTS ia32)
    endif()

    # Set default components
    set(MKL_FORTRAN_DEFAULT gf)
    set(MKL_THREAD_DEFAULT sequential)
    set(MKL_ARCH_DEFAULT lp64)
    set(MKL_MATH_DEFAULT blas lapack)


    # Get enabled languages
    get_property(LANG GLOBAL PROPERTY ENABLED_LANGUAGES)
    if(C IN_LIST LANG)
        set(C ${C})
    endif()
    if(CXX IN_LIST LANG)
        set(CXX ${CXX})
    endif()
    if(Fortran IN_LIST LANG)
        set(Fortran ${Fortran})
    endif()

    # Define required components
    foreach(fort ${MKL_FORTRAN_VARIANTS})
        if(${fort} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_FORTRAN_COMPONENT ${fort})
        endif()
    endforeach()

    foreach(arch ${MKL_ARCH_VARIANTS})
        if(${arch} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_ARCH_COMPONENT ${arch})
        endif()
    endforeach()

    foreach(thread ${MKL_THREAD_VARIANTS})
        if(${thread} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_THREAD_COMPONENT ${thread})
        endif()
    endforeach()
    foreach(math ${MKL_MATH_VARIANTS})
        if(${math} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_MATH_COMPONENTS ${math})
        endif()
    endforeach()

    #  If no component was requested, take the default
    if(NOT MKL_FIND_FORTRAN_COMPONENT)
        set(MKL_FIND_FORTRAN_COMPONENT ${MKL_FORTRAN_DEFAULT})
    endif()
    if(NOT MKL_FIND_THREAD_COMPONENT)
        set(MKL_FIND_THREAD_COMPONENT ${MKL_THREAD_DEFAULT})
    endif()
    if(NOT MKL_FIND_ARCH_COMPONENT)
        set(MKL_FIND_ARCH_COMPONENT ${MKL_ARCH_DEFAULT})
    endif()
    if(NOT MKL_FIND_MATH_COMPONENTS)
        set(MKL_FIND_MATH_COMPONENTS ${MKL_MATH_DEFAULT})
    endif()

    # If multiple components matched, the resulting target will not be unique, which means
    # the component list is ill-defined
    list(LENGTH MKL_FIND_FORTRAN_COMPONENT MKL_FIND_FORTRAN_LENGTH)
    list(LENGTH MKL_FIND_THREAD_COMPONENT MKL_FIND_THREAD_LENGTH)
    list(LENGTH MKL_FIND_ARCH_COMPONENT MKL_FIND_ARCH_LENGTH)
    if(MKL_FIND_FORTRAN_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Fortran components: ${MKL_FIND_FORTRAN_COMPONENT}\n")
    endif()
    if(MKL_FIND_THREAD_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Thread components: ${MKL_FIND_THREAD_COMPONENT}\n")
    endif()
    if(MKL_FIND_ARCH_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Architechture components: ${MKL_FIND_ARCH_COMPONENT}\n")
    endif()
    if(FAIL_MESSAGE)
        message(FATAL_ERROR "FindMKL component error: \n${FAIL_MESSAGE}")
    endif()


    # Define usable targets
    foreach (fort ${MKL_FIND_FORTRAN_COMPONENT})
        foreach (thread ${MKL_FIND_THREAD_COMPONENT})
            foreach (arch ${MKL_FIND_ARCH_COMPONENT})
                if(arch MATCHES "32" AND NOT TARGET mkl::mkl_${fort}_${arch})
                    # We make an alias for the mkl::mkl_gf library as mkl::mkl_gf_ia32
                    add_library(mkl::mkl_${fort}_${arch} ALIAS mkl::mkl_${fort})
                endif()
                add_library(mkl::mkl_${fort}_${thread}_${arch} INTERFACE IMPORTED)
                if(NOT MSVC)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE -Wl,--no-as-needed)
                endif()
                if(blas IN_LIST MKL_FIND_MATH_COMPONENTS)
                    target_link_libraries(mkl::mkl_blas${MKL_95_SUFFIX}${arch} INTERFACE -Wl,--start-group)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_blas${MKL_95_SUFFIX}${arch})
                    set(MKL_blas_FOUND TRUE PARENT_SCOPE)
                endif()

                if(lapack IN_LIST MKL_FIND_MATH_COMPONENTS)
                    target_link_libraries(mkl::mkl_lapack${MKL_95_SUFFIX}${arch} INTERFACE -Wl,--start-group)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_lapack${MKL_95_SUFFIX}${arch})
                    set(MKL_lapack_FOUND TRUE PARENT_SCOPE)
                endif()

                if(MSVC)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE
                            mkl::mkl_${fort}_${arch}
                            mkl::mkl_${thread}
                            mkl::mkl_core
                            )
                else()
                    target_link_libraries(mkl::mkl_${fort}_${arch} INTERFACE -Wl,--end-group)
                    target_link_libraries(mkl::mkl_${thread} INTERFACE -Wl,--end-group)
                    target_link_libraries(mkl::mkl_core INTERFACE -Wl,--end-group)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE
                            -Wl,--start-group
                            mkl::mkl_${fort}_${arch}
                            mkl::mkl_${thread}
                            mkl::mkl_core
                            -Wl,--end-group
                            -Wl,--as-needed
                            )
                endif()
                target_include_directories(mkl::mkl_${fort}_${thread}_${arch} SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
                set_target_properties(mkl::mkl_${fort}_${thread}_${arch} PROPERTIES INTERFACE_LINK_DIRECTORIES ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
                if(arch MATCHES "64")
                    target_compile_options(mkl::mkl_${fort}_${thread}_${arch} INTERFACE -m64)
                endif()
                if(thread MATCHES "gnu_thread")
                    find_package(Threads REQUIRED)
                    find_package(OpenMP COMPONENTS ${C} ${CXX} ${Fortran} REQUIRED)
                    foreach(lang ${LANG})
                        if(TARGET OpenMP::OpenMP_${lang})
                            target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE OpenMP::OpenMP_${lang})
                        endif()
                    endforeach()
                endif()
                if(thread MATCHES "intel_thread")
                    find_package(Threads REQUIRED)
                    target_link_libraries(mkl::iomp5 INTERFACE Threads::Threads)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::iomp5)
                endif()
                if(thread MATCHES "tbb_thread")
                    find_package(Threads REQUIRED)
                    target_link_libraries(mkl::tbb INTERFACE Threads::Threads)
                    target_link_libraries(mkl::tbbmalloc INTERFACE Threads::Threads)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::tbb mkl::tbbmalloc)
                endif()
                if(thread MATCHES "thread")
                    find_package(Threads REQUIRED)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE Threads::Threads)
                endif()
                target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE m dl)
                list(APPEND MKL_TARGETS mkl::mkl_${fort}_${thread}_${arch})
                set(MKL_${fort}_FOUND   TRUE PARENT_SCOPE)
                set(MKL_${thread}_FOUND TRUE PARENT_SCOPE)
                set(MKL_${arch}_FOUND   TRUE PARENT_SCOPE)
            endforeach()
        endforeach()
    endforeach()

    set(MKL_TARGETS ${MKL_TARGETS} PARENT_SCOPE)
endfunction()

# Test MKL
function(check_mkl_compiles)
    include(CheckCXXSourceCompiles)
    foreach(tgt ${MKL_TARGETS})
        unset(CMAKE_REQUIRED_LIBRARIES)
        if(NOT BUILD_SHARED_LIBS)
            set(CMAKE_REQUIRED_LIBRARIES -static-libgcc -static-libstdc++)
        endif()
        string(SUBSTRING ${tgt} 5 -1 tgt_name)
        list(APPEND CMAKE_REQUIRED_LIBRARIES ${tgt})
        if(TARGET ${tgt} AND ${tgt} MATCHES "sequential")
            check_cxx_source_compiles("
                        #include <mkl.h>
                        int main() {
                            const MKL_INT nx = 10, incx = 1, incy = 1;
                            double x[10], y[10];
                            for(int i = 0; i < 10; i++) x[i] = double(i);
                            dcopy(&nx, x, &incx, y, &incy);
                            return 0;
                        }
                        " COMPILES_${tgt_name})
        else()
            check_cxx_source_compiles("
                        #include <mkl.h>
                        int main() {
                            mkl_set_num_threads(2);
                            const MKL_INT nx = 10, incx = 1, incy = 1;
                            double x[10], y[10];
                            for(int i = 0; i < 10; i++) x[i] = double(i);
                            dcopy(&nx, x, &incx, y, &incy);
                            return 0;
                        }
                        " COMPILES_${tgt_name})
        endif()
        if(NOT COMPILES_${tgt_name})
            if(EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
                file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
                message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
            endif()
            message(FATAL_ERROR "Unable to compile a simple MKL program")
        endif()
    endforeach()
endfunction()



find_mkl_libraries()
setup_mkl_targets()
check_mkl_compiles()

list(LENGTH MKL_TARGETS MKL_TARGETS_NUM)
if(MKL_TARGETS_NUM EQUAL 1)
    add_library(mkl::mkl INTERFACE IMPORTED)
    target_link_libraries(mkl::mkl INTERFACE ${MKL_TARGETS})

    if(blas IN_LIST MKL_FIND_COMPONENTS)
        add_library(BLAS::BLAS INTERFACE IMPORTED)
        target_link_libraries(BLAS::BLAS INTERFACE ${MKL_TARGETS})
    endif()
    if(lapack IN_LIST MKL_FIND_COMPONENTS)
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
        target_link_libraries(LAPACK::LAPACK INTERFACE ${MKL_TARGETS})
    endif()
else()
    message(WARNING "Multiple MKL libraries found: ${MKL_TARGETS}"
                    "Try specifying unique COMPONENTS"
            )
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL
        REQUIRED_VARS MKL_INCLUDE_DIR MKL_TARGETS
        REASON_FAILURE_MESSAGE ${FAIL_MESSAGE}
        HANDLE_COMPONENTS
        )
