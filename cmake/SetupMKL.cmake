

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# or define -DEIGEN_USE_MKL_ALL.
# REMEMBER TO SET LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64 to use the shared libraries
# You can set that either as an environment variable in clion or add it to
# sudo nano /etc/ld.so.conf.d/intel_mkl.conf
# that simply contains the path to MKL libraries: /opt/intel/mkl/lib/intel64. Finally, call
# sudo ldconfig
# Acutally, sometimes all you need to do is sudo ldconfig to update the linker, because intel mkl installation
# automatically puts a file in /etc/ld.so.conf.d



#########################
# Getting SIGSEGVS / SEGFAULTS?
# When linking statically, pthreads need to be linked in a special way:
# set(PTHREAD_LIBRARY "-Wl,--whole-archive  -lpthread -Wl,--no-whole-archive") # This is needed with static linking!
# Use the quotations to force -lpthread to remain inside the whole-archive directive, otherwise cmake black magic
# can put it outside. Check with cmake verbose ON that the linking goes through as above.
#########################

#    set(MKL_USE_STATIC_LIBS ON)

if(NOT TARGET mkl::mkl AND DMRG_ENABLE_MKL)
    if(DMRG_ENABLE_OPENMP)
        find_package(OpenMP REQUIRED)
    endif()
    set(MKL_THREAD_DEFAULT none)
    if(TARGET openmp::openmp)
        if(OMP_LIBNAME MATCHES iomp)
            set(MKL_THREAD_DEFAULT intel)
        elseif(OMP_LIBNAME MATCHES gomp|libomp|CXX)
            set(MKL_THREAD_DEFAULT gnu)
        endif()
    endif()
    message(STATUS "MKL_THREAD_DEFAULT: ${MKL_THREAD_DEFAULT}")
    message(STATUS "OMP_LIBNAME       : ${OMP_LIBNAME}")
    # Make an "enum" for valid download methods
    set(MKL_THREADING_VALID none gnu intel tbb)
    set(MKL_THREADING_LAYER ${MKL_THREAD_DEFAULT} CACHE STRING "Intel MKL Threading Layer")
    set_property(CACHE MKL_THREADING_LAYER PROPERTY STRINGS ${MKL_THREADING_VALID})
    if (NOT MKL_THREADING_LAYER IN_LIST MKL_THREADING_VALID)
        message(FATAL_ERROR "MKL_THREADING_LAYER must be one of ${MKL_THREADING_VALID}")
    endif ()
    if(NOT BUILD_SHARED_LIBS AND MKL_THREADING_LAYER MATCHES tbb)
        message(FATAL_ERROR "MKL Threading library [tbb] requires dynamic linking. Use BUILD_SHARED_LIBS:BOOL=ON")
    endif()

    set(MKL_USE_SINGLE_DYNAMIC_LIBRARY OFF) # This doesn't work for some reason... You need to use the mkl_set_interface_layer(int) to select at runtime, which is not good when building dependencies!
    if (MKL_USE_SINGLE_DYNAMIC_LIBRARY AND NOT BUILD_SHARED_LIBS)
        message(WARNING "Disabling single dynamic mkl library\nCan't use MKL_USE_SINGLE_DYNAMIC_LIBRARY and -static simultaneously.")
    endif()

    find_package(MKL)
    if (NOT MKL_FOUND)
        message(WARNING "\
            Could not find Intel MKL library as requested by
            passing \"-DDMRG_ENABLE_MKL:BOOL=ON\" to CMake.")
    endif()



    if (MKL_FOUND)
        # Make a handle library for convenience. This "mkl" library is available throughout this cmake project later.
        add_library(mkl::mkl INTERFACE IMPORTED)
        target_compile_definitions(mkl::mkl INTERFACE MKL_AVAILABLE)
        if(MKL_THREADING_LAYER MATCHES gnu)
            find_package(OpenMP REQUIRED)
            find_package(Fortran REQUIRED)
            target_link_libraries(mkl::mkl INTERFACE mkl::mkl_gf_lp_gthread)
            get_target_property(MKL_INCLUDE_DIR mkl::mkl_gf_lp_gthread INTERFACE_INCLUDE_DIRECTORIES)
            target_include_directories(mkl::mkl SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
            target_link_libraries(mkl::mkl INTERFACE  gfortran::gfortran openmp::openmp)
        elseif(MKL_THREADING_LAYER MATCHES intel)
            find_package(Fortran REQUIRED)
            find_package(Threads REQUIRED)
            target_link_libraries(mkl::mkl INTERFACE mkl::mkl_gf_lp_ithread gfortran::gfortran mkl::iomp5 Threads::Threads m dl)
            get_target_property(MKL_INCLUDE_DIR mkl::mkl_gf_lp_ithread INTERFACE_INCLUDE_DIRECTORIES)
            target_include_directories(mkl::mkl SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
        elseif(MKL_THREADING_LAYER MATCHES tbb)
            find_package(Threads REQUIRED)
            target_link_libraries(mkl::mkl INTERFACE mkl::mkl_intel_lp_tthread mkl::tbb Threads::Threads m dl)
            get_target_property(MKL_INCLUDE_DIR mkl::mkl_intel_lp_tthread INTERFACE_INCLUDE_DIRECTORIES)
            target_include_directories(mkl::mkl SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
        elseif(MKL_THREADING_LAYER MATCHES none)
            find_package(Threads REQUIRED)
            find_package(Fortran REQUIRED)
            target_link_libraries(mkl::mkl INTERFACE mkl::mkl_gf_lp_seq)
            get_target_property(MKL_INCLUDE_DIR mkl::mkl_gf_lp_seq INTERFACE_INCLUDE_DIRECTORIES)
            target_include_directories(mkl::mkl SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
            target_link_libraries(mkl::mkl INTERFACE gfortran::gfortran Threads::Threads)
        endif()

        # Make the rest of the build structure aware of blas and lapack included in MKL.
        add_library(BLAS::BLAS          INTERFACE IMPORTED)
        add_library(LAPACK::LAPACK      INTERFACE IMPORTED)
        target_link_libraries(BLAS::BLAS        INTERFACE mkl::mkl)
        target_link_libraries(LAPACK::LAPACK    INTERFACE mkl::mkl)

        #   Test MKL
        function(check_mkl_compiles TARGETS)
            if(NOT BUILD_SHARED_LIBS)
                list(APPEND CMAKE_REQUIRED_LIBRARIES -static)
            endif()
            list(APPEND CMAKE_REQUIRED_LIBRARIES     ${TARGETS})
            include(CheckCXXSourceCompiles)
            check_cxx_source_compiles("
                #include <mkl.h>
                int main() {
                    const MKL_INT nx = 10, incx = 1, incy = 1;
                    double x[10], y[10];
                    for(int i = 0; i < 10; i++) x[i] = double(i);
                    dcopy(&nx, x, &incx, y, &incy);
                    return 0;
                }
                " MKL_COMPILES)

            if(NOT MKL_COMPILES)
                unset(MKL_COMPILES CACHE)
                unset(MKL_COMPILES PARENT_SCOPE)
                if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
                    file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
                    message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
                endif()
                message(FATAL_ERROR "Unable to compile a simple MKL program")
            endif()
                                #include <omp.h>

            if(NOT MKL_THREADING_LAYER MATCHES none)
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
                    " MKL_THREAD_COMPILES)
                if(NOT MKL_THREAD_COMPILES)
                    unset(MKL_THREAD_COMPILES CACHE)
                    unset(MKL_THREAD_COMPILES PARENT_SCOPE)
                    if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
                        file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
                        message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
                    endif()
                    message(FATAL_ERROR "Unable to compile a threaded MKL program")
                endif()
            endif()
        endfunction()
        check_mkl_compiles(mkl::mkl)
    endif()
endif()
