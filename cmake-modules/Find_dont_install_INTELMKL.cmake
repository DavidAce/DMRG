

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

if (USE_MKL)
    #    set(MKL_USE_STATIC_LIBS ON)
    set(MKL_MULTI_THREADED ${USE_OpenMP})
    set(MKL_USE_SINGLE_DYNAMIC_LIBRARY OFF) # This doesn't work for some reason... You need to use the mkl_set_interface_layer(int) to select at runtime, which is not good when building dependencies!
    if (MKL_USE_SINGLE_DYNAMIC_LIBRARY AND NOT BUILD_SHARED_LIBS)
        message(WARNING "Disabling single dynamic mkl library\nCan't use MKL_USE_SINGLE_DYNAMIC_LIBRARY and -static simultaneously.")
    endif()

    find_package(MKL)
    if (NOT MKL_FOUND)
        message(WARNING "\
        Could not find Intel MKL library. Turn off this
        warning by passing \"-DUSE_MKL:BOOL=OFF\" to CMake.")
    endif()

endif()



if (MKL_FOUND)


    #The order of these libraries is important when doing static linking!
    #To find out the order, check the Intel link line advisor.
    if(BUILD_SHARED_LIBS AND MKL_USE_SINGLE_DYNAMIC_LIBRARY)
        #This doesn't seem to work,  The gnu fortran GF_LP library doesn't get enabled --> arpack++ complains!
        set(MKL_LIBRARIES -Wl,--no-as-needed ${MKL_RT_LIBRARY} -Wl,--as-needed  )
    else()
        set(MKL_LIBRARIES  -Wl,--no-as-needed ${MKL_BLAS_LP_LIBRARY} ${MKL_LAPACK_LP_LIBRARY}   -Wl,--start-group  ${MKL_GF_LP_LIBRARY})

        if(MKL_MULTI_THREADED)
            list(APPEND MKL_LIBRARIES  ${MKL_GNUTHREAD_LIBRARY} ${MKL_INTELTHREAD_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group )
            if(BUILD_SHARED_LIBS)
                list(APPEND MKL_LIBRARIES ${MKL_IOMP5_LIBRARY})
            else()
                list(APPEND MKL_LIBRARIES ${OpenMP_LIBRARIES})
            endif()
        else()
            list(APPEND MKL_LIBRARIES  ${MKL_SEQUENTIAL_LIBRARY}  ${MKL_CORE_LIBRARY}  -Wl,--end-group -Wl,--as-needed )
        endif()
    endif()

    set(MKL_FLAGS -m64 -I${MKL_ROOT_DIR}/lib/intel64/lp64 )
    if(BUILD_SHARED_LIBS)
        list(APPEND MKL_FLAGS -fPIC)
    endif()

    add_definitions(-DMKL_AVAILABLE)
    #    add_definitions(-DMKL_VERBOSE -DMKL_Complex8=std::complex<float> -DMKL_Complex16=std::complex<double>)

    # Make a handle library for convenience. This "mkl" library is available throughout this cmake project later.
    add_library(mkl INTERFACE)
    list(APPEND MKL_LIBRARIES Threads::Threads)
    target_link_libraries(mkl INTERFACE ${MKL_LIBRARIES}  -ldl -lm gfortran)
    target_include_directories(mkl INTERFACE ${MKL_INCLUDE_DIR})
    target_compile_options(mkl INTERFACE ${MKL_FLAGS})
    set_target_properties(mkl PROPERTIES INTERFACE_LINK_DIRECTORIES  "${MKL_ROOT_DIR}/lib/intel64")
    # BLAS and LAPACK are included in the MKL.
    set(BLAS_LIBRARIES   ${MKL_LIBRARIES})
    set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
    set(FC_LDLAGS -lm -ldl ${GFORTRAN_LIB}) #This one is needed if any sub projects wants to link its own stuff using MKL. For instance, arpack-ng.
    #    set(FC_LDLAGS -lm -ldl  -fPIC ${PTHREAD_LIBRARY}) #This one is needed if any sub projects wants to link its own stuff using MKL. For instance, arpack-ng.


    # Make the rest of the build structure aware of blas and lapack included in MKL.
    add_library(blas ALIAS mkl)
    add_library(lapack ALIAS mkl)
    #    add_library(lapacke ALIAS mkl)

    message("")
    message("============================ MKL SUMMARY ===================================")
    message("MKL_LIBRARIES                  : ${MKL_LIBRARIES}" )
    message("MKL_RT_LIBRARY                 : ${MKL_RT_LIBRARY}" )
    message("MKL_INCLUDE_DIR                : ${MKL_INCLUDE_DIR}" )
    message("MKL_FLAGS                      : ${MKL_FLAGS}" )
    message("MKLROOT                        : $ENV{MKLROOT}" )
    message("MKL_USE_SINGLE_DYNAMIC_LIBRARY : ${MKL_USE_SINGLE_DYNAMIC_LIBRARY}" )
    message("=============================================================================")
    message("")


    #   Test features
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES ${MKL_LIBRARIES}  ${FC_LDLAGS})
    set(CMAKE_REQUIRED_INCLUDES  ${MKL_INCLUDE_DIR})
    set(CMAKE_REQUIRED_FLAGS     ${MKL_FLAGS})
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
        message(FATAL_ERROR "Unable to compile a simple MKL program")
    endif()

    if(MKL_MULTI_THREADED)
        check_cxx_source_compiles("
            #include <omp.h>
            #include <mkl.h>
            int main() {
                mkl_set_num_threads(2);
                const MKL_INT nx = 10, incx = 1, incy = 1;
                double x[10], y[10];
                for(int i = 0; i < 10; i++) x[i] = double(i);
                dcopy(&nx, x, &incx, y, &incy);
                return 0;
            }
            " MKL_OMP_COMPILES)

        if(NOT MKL_OMP_COMPILES)
            message(FATAL_ERROR "Unable to compile a simple MKL program with OMP")
        endif()
    endif()

endif()


#
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    if("${_variableName}" MATCHES "BLAS"
#            OR "${_variableName}" MATCHES "blas"
#            OR "${_variableName}" MATCHES "LAPACK"
#            OR "${_variableName}" MATCHES "lapack")
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endif()
#endforeach()

#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        if("${_variableName}" MATCHES "MKL" OR "${_variableName}" MATCHES "mkl")
#            message(STATUS "${_variableName}=${${_variableName}}")
#        endif()
#    endforeach()
