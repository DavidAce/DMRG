
function(strip_genex input_string output_string)
    set(result_string)
    foreach(elem ${input_string})
        if(${elem} MATCHES "[$]")
            string(REGEX MATCHALL "([-=]+[a-z0-9]+)" filtered_string "${input_string}")
            list(APPEND result_string ${filtered_string})
        else()
            list(APPEND result_string ${elem})
        endif()
    endforeach()
    set(${output_string} ${result_string} PARENT_SCOPE)
endfunction()


function(check_omp_compiles omp_tgt)
    if(NOT TARGET ${omp_tgt})
        message(FATAL_ERROR "Given OpenMP target [${omp_tgt}] is not a target")
    endif()
    include(CheckIncludeFileCXX)
    include(cmake-modules/getExpandedTarget.cmake)
    expand_target_libs("${omp_tgt}" omp_lib)
    expand_target_incs("${omp_tgt}" omp_inc)
    expand_target_opts("${omp_tgt}" omp_opt)
    expand_target_defs("${omp_tgt}" omp_def)
    if(NOT BUILD_SHARED_LIBS)
        list(APPEND CMAKE_REQUIRED_LIBRARIES -static)
    endif()
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${omp_tgt}) # Can be a ;list
    if(DMRG_PRINT_CHECKS)
        message(STATUS "OPENMP COMPILE TEST required      ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "OPENMP COMPILE TEST target        ${omp_tgt}")
        message(STATUS "OPENMP COMPILE TEST libraries     ${omp_lib}")
        message(STATUS "OPENMP COMPILE TEST options       ${omp_opt}")
        message(STATUS "OPENMP COMPILE TEST includes      ${omp_inc}")
        message(STATUS "OPENMP COMPILE TEST definitions   ${omp_def}")
    endif()
    unset(has_omp_h)
    unset(has_omp_h CACHE)
    unset(OMP_COMPILES)
    unset(OMP_COMPILES CACHE)

    check_include_file_cxx(omp.h    has_omp_h)
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
            #include <omp.h>
            #include <iostream>
            #ifndef _OPENMP
            #error You forgot to define _OPENMP
            #endif
            int main() {
                omp_set_num_threads(4);
                std::cout << \"OMP Threads \" << omp_get_max_threads() << std::endl;
                int counter = 0;
                #pragma omp parallel for shared(counter)
                for(int i = 0; i < 16; i++){
                    #pragma omp atomic
                    counter++;
                }
                std::cout << \"Counter is: \" << counter << std::endl;

                return 0;
            }
            " OMP_COMPILES
            )
    set(OMP_COMPILES ${OMP_COMPILES} PARENT_SCOPE)
    if(NOT OMP_COMPILES)
        unset(OMP_COMPILES CACHE)
        unset(OMP_COMPILES PARENT_SCOPE)
    endif()

endfunction()





if(NOT OpenMP_FOUND AND NOT TARGET openmp::openmp AND BUILD_SHARED_LIBS)
    set(BACKUP ${CMAKE_MODULE_PATH})
    unset(CMAKE_MODULE_PATH)
    find_package(OpenMP)
    list(APPEND CMAKE_MODULE_PATH ${BACKUP})
    unset(BACKUP)
    if(OpenMP_FOUND AND TARGET OpenMP::OpenMP_CXX)
        set_target_properties(OpenMP::OpenMP_CXX PROPERTIES INTERFACE_COMPILE_DEFINITIONS "_OPENMP")
        check_omp_compiles(OpenMP::OpenMP_CXX)
        if(OMP_COMPILES AND TARGET OpenMP::OpenMP_CXX)
            add_library(openmp::openmp INTERFACE IMPORTED)
            target_link_libraries(openmp::openmp INTERFACE OpenMP::OpenMP_CXX)
        endif()
    endif()
endif()



if(NOT OpenMP_FOUND)
    find_package(Threads)
    list(APPEND omp_lib_candidates gomp omp iomp5)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # MKL comes with a useable iomp5 library
        find_path(MKL_ROOT_DIR
                include/mkl.h
                HINTS ${CMAKE_INSTALL_PREFIX}
                PATHS
                $ENV{MKL_DIR}  ${MKL_DIR}
                $ENV{MKLDIR}   ${MKLDIR}
                $ENV{MKLROOT}  ${MKLROOT}
                $ENV{MKL_ROOT} ${MKL_ROOT}
                $ENV{mkl_root} ${mkl_root}
                $ENV{HOME}/intel/mkl
                /opt/intel/mkl
                /opt/intel
                $ENV{BLAS_DIR}
                /usr/lib/x86_64-linux-gnu
                /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
                "Program Files (x86)/Intel/ComposerXE-2011/mkl"
                PATH_SUFFIXES
                intel intel/mkl mkl lib
                )

        find_library(omp_lib_iomp5 NAMES iomp5
                PATHS
                /usr/lib/llvm-11
                /usr/lib/llvm-10
                /usr/lib/llvm-9
                /usr/lib/llvm-8
                /usr/lib/llvm-7
                /usr
                ${MKL_ROOT_DIR}
                /opt/intel/lib/intel64
                PATH_SUFFIXES
                include lib local ../intel/lib/intel64
                )
        if(omp_lib_iomp5)
            list(APPEND omp_lib_candidates ${omp_lib_iomp5})
        endif()

    endif()


    foreach(lib ${omp_lib_candidates})
        # Make a target wrapper for a library
        if(NOT TARGET openmp::openmp)
            add_library(openmp::openmp INTERFACE IMPORTED)
        endif()
        if(IS_ABSOLUTE lib)
            set_target_properties(openmp::openmp PROPERTIES IMPORTED_LOCATION ${lib})
            set_target_properties(openmp::openmp PROPERTIES INTERFACE_LINK_LIBRARIES "Threads::Threads;rt;dl")
        else()
            set_target_properties(openmp::openmp PROPERTIES INTERFACE_LINK_LIBRARIES "${lib};Threads::Threads;rt;dl")
        endif()
        set_target_properties(openmp::openmp PROPERTIES INTERFACE_COMPILE_DEFINITIONS "_OPENMP")
        check_omp_compiles(openmp::openmp)
        if(OMP_COMPILES)
            break()
        endif()
    endforeach()
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMP DEFAULT_MSG OMP_COMPILES)


