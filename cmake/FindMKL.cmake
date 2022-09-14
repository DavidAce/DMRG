# find_package module for the Intel Math Kernel Library (MKL)
#
# COMPONENTS
# architecture:

# Once done this will define
#
# MKL_FOUND - system has MKL
# MKL_ROOT_DIR - path to the MKL base directory
# MKL_INCLUDE_DIR - the MKL include directory
# MKL_TARGETS - MKL CMake targets
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

    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES .so.1 .so.2) # For tbb

    # macos
    if(APPLE)
        set(MKL_ARCH_DIR "intel64")
    endif()

    if(FORCE_BUILD_32BITS)
        set(MKL_ARCH_DIR "ia32")
    endif()

    if(WIN32)
        if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
            set(MKL_ARCH_DIR "intel64")
        else()
            set(MKL_ARCH_DIR "ia32")
        endif()
    endif()
    set(MKL_ARCH_DIR ${MKL_ARCH_DIR} PARENT_SCOPE)

    set(MKL_ROOT_SEARCH_PATHS
        ${MKL_ROOT_DIR}
        $ENV{MKL_DIR} ${MKL_DIR}
        $ENV{MKLDIR} ${MKLDIR}
        $ENV{MKLROOT} ${MKLROOT}
        $ENV{MKL_ROOT} ${MKL_ROOT}
        $ENV{mkl_root} ${mkl_root}
        $ENV{HOME}
        /opt
        /opt/intel
        /opt/intel/oneapi
        /usr/lib/x86_64-linux-gnu
        /usr
        /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
        "Program Files (x86)/Intel/ComposerXE-2011/mkl"
        )

    set(MKL_PATH_SUFFIXES
        intel/oneapi/mkl/latest
        oneapi/mkl/latest
        intel/oneapi/mkl
        oneapi/mkl
        intel/mkl
        mkl
        )

    find_path(MKL_ROOT_DIR
              include/mkl.h
              PATHS ${MKL_ROOT_SEARCH_PATHS}
              PATH_SUFFIXES ${MKL_PATH_SUFFIXES}
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
        set(mkl_libnames rt core sequential intel_thread gnu_thread tbb_thread cdft_core)
        if(MKL_ARCH_DIR MATCHES "32")
            list(APPEND mkl_libnames
                 intel
                 gf
                 blas
                 lapack
                 )
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
                 scalapack_lp64
                 scalapack_ilp64
                 blacs_intelmpi_lp64
                 blacs_intelmpi_ilp64
                 blacs_openmpi_lp64
                 blacs_openmpi_ilp64
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
                message(TRACE "Adding target mkl::mkl_${lib} : ${MKL_${lib}_LIBRARY}")
            else()
                message(TRACE "Missed target mkl::mkl_${lib}")

            endif()
        endforeach()
        if(TARGET mkl::mkl_rt)
            set_target_properties(mkl::mkl_rt PROPERTIES INTERFACE_LINK_DIRECTORIES ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR})
        endif()

        foreach(lib ${par_libnames})
            find_library(MKL_${lib}_LIBRARY
                         ${lib}
                         HINTS
                         ${MKL_ROOT_DIR}
                         ${MKL_ROOT_DIR}/tbb/lib/${MKL_ARCH_DIR}
                         ${MKL_ROOT_DIR}/../tbb/lib/${MKL_ARCH_DIR}
                         ${MKL_ROOT_DIR}/../../tbb/lib/${MKL_ARCH_DIR}
                         ${MKL_ROOT_DIR}/tbb/latest/lib/${MKL_ARCH_DIR}
                         ${MKL_ROOT_DIR}/../tbb/latest/lib/${MKL_ARCH_DIR}
                         ${MKL_ROOT_DIR}/../../tbb/latest/lib/${MKL_ARCH_DIR}
                         PATH_SUFFIXES
                         lib lib/${MKL_ARCH_DIR}
                         gcc4.4
                         gcc4.7
                         gcc4.8
                         NO_DEFAULT_PATH
                         )
            if(MKL_${lib}_LIBRARY)
                add_library(mkl::${lib} UNKNOWN IMPORTED)
                set_target_properties(mkl::${lib} PROPERTIES IMPORTED_LOCATION "${MKL_${lib}_LIBRARY}")
                set_target_properties(mkl::${lib} PROPERTIES LINK_WHAT_YOU_USE TRUE)
                message(TRACE "Adding library mkl::${lib} : ${MKL_${lib}_LIBRARY}")
            else()
                message(TRACE "Missed library mkl::${lib}")
            endif()
        endforeach()
    endif()
endfunction()

function(setup_mkl_targets)

    # Setup all the library variants
    set(MKL_FORTRAN_VARIANTS intel gf)
    set(MKL_THREAD_VARIANTS sequential intel_thread gnu_thread tbb_thread)
    set(MKL_ARCH_VARIANTS)
    set(MKL_MATH_VARIANTS blas lapack)
    set(MKL_BLACS_VARIANTS blacs_intelmpi blacs_openmpi)
    set(MKL_CLUSTER_LIBS scalapack cdft)
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
    set(MKL_BLACS_DEFAULT)
    set(MKL_CLUSTER_DEFAULT)

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
            list(APPEND MKL_FIND_FORTRAN_COMPONENTS ${fort})
        endif()
    endforeach()

    foreach(arch ${MKL_ARCH_VARIANTS})
        if(${arch} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_ARCH_COMPONENTS ${arch})
        endif()
    endforeach()

    foreach(thread ${MKL_THREAD_VARIANTS})
        if(${thread} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_THREAD_COMPONENTS ${thread})
        endif()
    endforeach()
    foreach(math ${MKL_MATH_VARIANTS})
        if(${math} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_MATH_COMPONENTS ${math})
        endif()
    endforeach()

    foreach(blacs ${MKL_BLACS_VARIANTS})
        if(${blacs} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_BLACS_COMPONENTS ${blacs})
        endif()
    endforeach()

    foreach(libs ${MKL_CLUSTER_LIBS})
        if(${libs} IN_LIST MKL_FIND_COMPONENTS)
            list(APPEND MKL_FIND_CLUSTER_COMPONENTS ${libs})
        endif()
    endforeach()

    #  If no component was requested, take the default
    if(NOT MKL_FIND_FORTRAN_COMPONENTS)
        set(MKL_FIND_FORTRAN_COMPONENTS ${MKL_FORTRAN_DEFAULT})
    endif()
    if(NOT MKL_FIND_THREAD_COMPONENTS)
        set(MKL_FIND_THREAD_COMPONENTS ${MKL_THREAD_DEFAULT})
    endif()
    if(NOT MKL_FIND_ARCH_COMPONENTS)
        set(MKL_FIND_ARCH_COMPONENTS ${MKL_ARCH_DEFAULT})
    endif()
    if(NOT MKL_FIND_MATH_COMPONENTS)
        set(MKL_FIND_MATH_COMPONENTS ${MKL_MATH_DEFAULT})
    endif()
    if(NOT MKL_FIND_BLACS_COMPONENTS)
        set(MKL_FIND_BLACS_COMPONENTS ${MKL_BLACS_DEFAULT})
    endif()
    if(NOT MKL_FIND_CLUSTER_COMPONENTS)
        set(MKL_FIND_CLUSTER_COMPONENTS ${MKL_CLUSTER_DEFAULT})
    endif()

    # If multiple components matched, the resulting target will not be unique, which means
    # the component list is ill-defined
    list(LENGTH MKL_FIND_FORTRAN_COMPONENTS MKL_FIND_FORTRAN_LENGTH)
    list(LENGTH MKL_FIND_THREAD_COMPONENTS MKL_FIND_THREAD_LENGTH)
    list(LENGTH MKL_FIND_ARCH_COMPONENTS MKL_FIND_ARCH_LENGTH)
    list(LENGTH MKL_FIND_BLACS_COMPONENTS MKL_FIND_BLACS_LENGTH)
    list(LENGTH MKL_FIND_CLUSTER_COMPONENTS MKL_FIND_CLUSTER_LENGTH)
    if(MKL_FIND_FORTRAN_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Fortran components: ${MKL_FIND_FORTRAN_COMPONENTS}\n")
    endif()
    if(MKL_FIND_THREAD_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Thread components: ${MKL_FIND_THREAD_COMPONENTS}\n")
    endif()
    if(MKL_FIND_ARCH_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple Architechture components: ${MKL_FIND_ARCH_COMPONENTS}\n")
    endif()
    if(MKL_FIND_BLACS_LENGTH GREATER 1)
        list(APPEND FAIL_MESSAGE "Matched multiple BLACS components: ${MKL_FIND_BLACS_COMPONENTS}\n")
    endif()
    if(MKL_FIND_CLUSTER_LENGTH GREATER 2)
        list(APPEND FAIL_MESSAGE "Matched too many Cluster components: ${MKL_FIND_CLUSTER_COMPONENTS}\n")
    endif()
    if(MKL_FIND_CLUSTER_LENGTH GREATER 0 AND MKL_FIND_BLACS_LENGTH EQUAL 0)
        list(APPEND FAIL_MESSAGE "Components [${MKL_FIND_CLUSTER_COMPONENTS}] requires choosing [${MKL_BLACS_VARIANTS}]\n")
    endif()
    if(FAIL_MESSAGE)
        message(FATAL_ERROR "FindMKL component error: \n${FAIL_MESSAGE}")
    endif()

    # Define usable targets
    foreach(fort ${MKL_FIND_FORTRAN_COMPONENTS})
        foreach(thread ${MKL_FIND_THREAD_COMPONENTS})
            foreach(arch ${MKL_FIND_ARCH_COMPONENTS})
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

                if(cdft IN_LIST MKL_FIND_CLUSTER_COMPONENTS)
                    target_link_libraries(mkl::mkl_cdft_${arch} INTERFACE -Wl,--end-group MPI::MPI_C)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_cdft_${arch})
                    set(MKL_cdft_FOUND TRUE PARENT_SCOPE)
                endif()

                if(scalapack IN_LIST MKL_FIND_CLUSTER_COMPONENTS)
                    target_link_libraries(mkl::mkl_scalapack_${arch} INTERFACE -Wl,--start-group MPI::MPI_C)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_scalapack_${arch})
                    set(MKL_scalapack_FOUND TRUE PARENT_SCOPE)
                endif()

                if(blacs_openmpi IN_LIST MKL_FIND_BLACS_COMPONENTS)
                    target_link_libraries(mkl::mkl_blacs_openmpi_${arch} INTERFACE -Wl,--end-group MPI::MPI_C)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_blacs_openmpi_${arch})
                    target_link_libraries(mkl::mkl_core INTERFACE mkl::mkl_blacs_openmpi_${arch})
                    if(TARGET mkl::mkl_scalapack_${arch})
                        target_link_libraries(mkl::mkl_scalapack_${arch} INTERFACE mkl::mkl_blacs_openmpi_${arch})
                    endif()
                    set(MKL_blacs_openmpi_FOUND TRUE PARENT_SCOPE)
                endif()

                if(blacs_intelmpi IN_LIST MKL_FIND_BLACS_COMPONENTS)
                    target_link_libraries(mkl::mkl_blacs_intelmpi_${arch} INTERFACE -Wl,--end-group MPI::MPI_C)
                    target_link_libraries(mkl::mkl_${fort}_${thread}_${arch} INTERFACE mkl::mkl_blacs_intelmpi_${arch})
                    target_link_libraries(mkl::mkl_core INTERFACE mkl::mkl_blacs_intelmpi_${arch})
                    if(TARGET mkl::mkl_scalapack_${arch})
                        target_link_libraries(mmkl::mkl_scalapack_${arch} INTERFACE mkl::mkl_blacs_intelmpi_${arch})
                    endif()
                    set(MKL_blacs_intelmpi_FOUND TRUE PARENT_SCOPE)
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
                    find_package(OpenMP COMPONENTS ${C} ${CXX} ${Fortran} REQUIRED)
                    foreach(lang ${LANG})
                        if(OpenMP_${lang}_FLAGS)
                            target_link_options(mkl::mkl_${fort}_${thread}_${arch} INTERFACE ${OpenMP_${lang}_FLAGS})
                            target_compile_options(mkl::mkl_${fort}_${thread}_${arch} INTERFACE
                                                   $<$<COMPILE_LANGUAGE:${LANG}>:${OpenMP_${lang}_FLAGS}>)
                        elseif(TARGET OpenMP::OpenMP_${lang})
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
                set(MKL_${fort}_FOUND TRUE PARENT_SCOPE)
                set(MKL_${thread}_FOUND TRUE PARENT_SCOPE)
                set(MKL_${arch}_FOUND TRUE PARENT_SCOPE)
            endforeach()
        endforeach()
    endforeach()

    set(MKL_TARGETS ${MKL_TARGETS} PARENT_SCOPE)
endfunction()

# Test MKL
function(check_mkl_compiles)
    include(CheckCXXSourceCompiles)
    foreach(tgt ${MKL_TARGETS})
        string(SUBSTRING ${tgt} 5 -1 tgt_name)
        set(CMAKE_REQUIRED_LIBRARIES ${tgt})
        include(cmake/PrintTargetInfo.cmake)
        print_target_info(${tgt} "")
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
            if(TB_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
                file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
                message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
            endif()
            message(FATAL_ERROR "Unable to compile a simple MKL program using target [${tgt}]. See error in ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
        endif()
    endforeach()
endfunction()

function(get_mkl_info libs incs opts defs ftrs target)
    if(NOT TARGET ${target})
        message(FATAL_ERROR "Given argument is not a valid CMake target: [${target}]")
    endif()
    get_target_property(LIB ${target} INTERFACE_LINK_LIBRARIES)
    get_target_property(INC ${target} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(OPT ${target} INTERFACE_COMPILE_OPTIONS)
    get_target_property(DEF ${target} INTERFACE_COMPILE_DEFINITIONS)
    get_target_property(FTR ${target} INTERFACE_COMPILE_FEATURES)
    get_target_property(TYP ${target} TYPE)
    get_target_property(IMP ${target} IMPORTED)
    if(IMP)
        if(NOT TYP MATCHES "INTERFACE" OR CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
            get_target_property(LOC ${target} LOCATION)
        endif()
    endif()
    if(LOC AND NOT TARGET ${LOC})
        list(APPEND ${libs} ${LOC})
    endif()
    if(INC)
        list(APPEND ${incs} ${INC})
    endif()
    if(OPT)
        list(APPEND ${opts} ${OPT})
    endif()
    if(DEF)
        list(APPEND ${defs} ${DEF})
    endif()
    if(FTR)
        list(APPEND ${ftrs} ${FTR})
    endif()
    foreach(lib ${LIB})
        if(TARGET ${lib})
            get_mkl_info(${libs} ${incs} ${opts} ${defs} ${ftrs} ${lib})
        else()
            list(APPEND ${libs} ${lib})
        endif()
    endforeach()

    list(REVERSE ${libs})
    list(REVERSE ${incs})
    list(REVERSE ${opts})
    list(REVERSE ${defs})
    list(REVERSE ${ftrs})

    list(REMOVE_DUPLICATES ${libs})
    list(REMOVE_DUPLICATES ${incs})
    list(REMOVE_DUPLICATES ${opts})
    list(REMOVE_DUPLICATES ${defs})
    list(REMOVE_DUPLICATES ${ftrs})

    list(REVERSE ${libs})
    list(REVERSE ${incs})
    list(REVERSE ${opts})
    list(REVERSE ${defs})
    list(REVERSE ${ftrs})

    set(${libs} ${${libs}} PARENT_SCOPE)
    set(${incs} ${${incs}} PARENT_SCOPE)
    set(${opts} ${${opts}} PARENT_SCOPE)
    set(${defs} ${${defs}} PARENT_SCOPE)
    set(${ftrs} ${${ftrs}} PARENT_SCOPE)

endfunction()

if(MKL_FOUND)
    return()
endif()

if(MKL_FIND_COMPONENTS MATCHES blacs|scalapack|cdft)
    find_package(MPI COMPONENTS C REQUIRED)
endif()

find_mkl_libraries()
setup_mkl_targets()
check_mkl_compiles()

list(LENGTH MKL_TARGETS MKL_TARGETS_NUM)
if(MKL_TARGETS_NUM GREATER 1)
    message(FATAL_ERROR "Multiple MKL libraries found: ${MKL_TARGETS}"
            "Try specifying unique COMPONENTS"
            )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL
                                  REQUIRED_VARS MKL_INCLUDE_DIR MKL_TARGETS
                                  REASON_FAILURE_MESSAGE ${FAIL_MESSAGE}
                                  HANDLE_COMPONENTS
                                  )

if(MKL_FOUND)
    if(NOT TARGET mkl::mkl)
        add_library(mkl::mkl INTERFACE IMPORTED)
        target_link_libraries(mkl::mkl INTERFACE ${MKL_TARGETS})
        if(blas IN_LIST MKL_FIND_COMPONENTS AND NOT TARGET BLAS::BLAS)
            add_library(BLAS::BLAS INTERFACE IMPORTED)
            target_link_libraries(BLAS::BLAS INTERFACE ${MKL_TARGETS})
        endif()
        if(lapack IN_LIST MKL_FIND_COMPONENTS AND NOT TARGET LAPACK::LAPACK)
            add_library(LAPACK::LAPACK INTERFACE IMPORTED)
            target_link_libraries(LAPACK::LAPACK INTERFACE ${MKL_TARGETS})
        endif()
    endif()
    unset(MKL_LIBRARIES)
    unset(MKL_LIBRARIES CACHE)
    get_mkl_info(MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_COMPILE_OPTIONS MKL_COMPILE_DEFINITIONS MKL_COMPILE_FEATURES ${MKL_TARGETS})

    set(MKL_LIBRARIES "${MKL_LIBRARIES}" CACHE INTERNAL "")
    set(MKL_INCLUDE_DIRS "${MKL_INCLUDE_DIRS}" CACHE INTERNAL "")
    set(MKL_COMPILE_OPTIONS "${MKL_COMPILE_OPTIONS}" CACHE INTERNAL "")
    set(MKL_COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}" CACHE INTERNAL "")
    set(MKL_COMPILE_FEATURES "${MKL_COMPILE_FEATURES}" CACHE INTERNAL "")

    # Do the same for blas
    set(BLAS_FOUND "TRUE" CACHE INTERNAL "")
    set(BLAS_LIBRARIES "${MKL_LIBRARIES}" CACHE INTERNAL "")
    set(BLAS_INCLUDE_DIRS "${MKL_INCLUDE_DIRS}" CACHE INTERNAL "")
    set(BLAS_COMPILE_OPTIONS "${MKL_COMPILE_OPTIONS}" CACHE INTERNAL "")
    set(BLAS_COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}" CACHE INTERNAL "")
    set(BLAS_COMPILE_FEATURES "${MKL_COMPILE_FEATURES}" CACHE INTERNAL "")

    # Do the same for lapack
    set(LAPACK_FOUND "TRUE" CACHE INTERNAL "")
    set(LAPACK_LIBRARIES "${MKL_LIBRARIES}" CACHE INTERNAL "")
    set(LAPACK_INCLUDE_DIRS "${MKL_INCLUDE_DIRS}" CACHE INTERNAL "")
    set(LAPACK_COMPILE_OPTIONS "${MKL_COMPILE_OPTIONS}" CACHE INTERNAL "")
    set(LAPACK_COMPILE_DEFINITIONS "${MKL_COMPILE_DEFINITIONS}" CACHE INTERNAL "")
    set(LAPACK_COMPILE_FEATURES "${MKL_COMPILE_FEATURES}" CACHE INTERNAL "")
    mark_as_advanced(MKL_LIBRARIES
                     MKL_INCLUDE_DIRS
                     MKL_COMPILE_OPTIONS
                     MKL_COMPILE_DEFINITIONS
                     MKL_COMPILE_FEATURES)

endif()