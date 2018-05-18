

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# Adding a -DEIGEN_USE_MKL_ALL here may conflict with arpack++

set(MKL_USE_STATIC_LIBS ON)
set(MKL_MULTI_THREADED OFF)
set(MKL_USE_SINGLE_DYNAMIC_LIBRARY OFF)
find_package(MKL)
if (MKL_FOUND)

    # Remove the libmkl_intel_lp64 library, which is only used when compiling with Intel Fortran.
    list(REMOVE_ITEM MKL_LIBRARIES ${MKL_INTEL_LP64_LIBRARY})
    list(INSERT MKL_LIBRARIES 0 ${MKL_GF_LP64_LIBRARY})
    list(REMOVE_DUPLICATES MKL_LIBRARIES)

    add_definitions(-DMKL_AVAILABLE)
    set(MKL_FLAGS -m64 -I${MKL_ROOT}/include)
    set(MKL_LFLAGS -L${MKL_ROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl)

    # Make a handle library for convenience. This "mkl" library is available throughout this cmake project later.
    add_library(mkl INTERFACE IMPORTED)
    set_target_properties(mkl PROPERTIES
            INTERFACE_LINK_LIBRARY "${MKL_LIBRARIES};${MKL_LFLAGS}"
            INTERFACE_INCLUDE_DIRECTORY "${MKL_INCLUDE_DIR}"
            INTERFACE_COMPILE_OPTIONS "${MKL_FLAGS}"
            )
    target_link_libraries(${PROJECT_NAME} PUBLIC mkl ${MKL_LIBRARIES} ${MKL_LFLAGS})
    target_compile_options(${PROJECT_NAME} PUBLIC ${MKL_FLAGS})
    target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INCLUDE_DIR})


    # BLAS and LAPACK are included in the MKL.
    set(BLAS_LIBRARIES ${MKL_LIBRARIES})
    set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
    set(EXTRA_LDLAGS ${MKL_LFLAGS}) #This one is needed if any sub projects wants to link its own stuff using MKL. For instance, arpack-ng.


    # Make the rest of the build structure aware of blas and lapack included in MKL.
    add_library(blas INTERFACE IMPORTED)
    set_target_properties(blas PROPERTIES
            INTERFACE_LINK_LIBRARIES        "${BLAS_LIBRARIES};${MKL_LFLAGS}"
            INTERFACE_INCLUDE_DIRECTORY     "${MKL_INCLUDE_DIR}"
            INTERFACE_COMPILE_OPTIONS       "${MKL_FLAGS}"
            )

    add_library(lapack INTERFACE IMPORTED)
    set_target_properties(lapack PROPERTIES
            INTERFACE_LINK_LIBRARIES        "${LAPACK_LIBRARIES};${MKL_LFLAGS}"
            INTERFACE_INCLUDE_DIRECTORY     "${MKL_INCLUDE_DIR}"
            INTERFACE_COMPILE_OPTIONS       "${MKL_FLAGS}"
            )


    #   Test features
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${MKL_INCLUDE_DIR})
    set(CMAKE_REQUIRED_LIBRARIES ${MKL_LIBRARIES} ${MKL_LFLAGS})
    set(CMAKE_REQUIRED_FLAGS ${MKL_FLAGS})
    check_cxx_source_compiles("
        #include <mkl.h>
        int main() {
            int nx = 10, incx = 1, incy = 1;
            double x[10], y[10];
            for(int i = 0; i < 10; i++) x[i] = double(i);
            dcopy(&nx, x, &incx, y, &incy);
            return 0;
        }
        " MKL_COMPILES)
    if(NOT MKL_COMPILES)
        message(FATAL_ERROR "Unable to compile a simple MKL program")
    endif(NOT MKL_COMPILES)

    message("")
    message("======MKL SUMMARY ======")
    message("MKL_LIBRARIES                             : ${MKL_LIBRARIES}" )
    message("MKL_RT_LIBRARY                            : ${MKL_RT_LIBRARY}" )
    message("MKL_INCLUDE_DIR                           : ${MKL_INCLUDE_DIR}" )
    message("MKL_FLAGS                                 : ${MKL_FLAGS}" )
    message("MKL_LFLAGS                                : ${MKL_LFLAGS}" )
    message("MKL_ROOT                                  : ${MKL_ROOT}" )
    message("MKLROOT                                   : ${MKLROOT}" )
    message("MKL_USE_SINGLE_DYNAMIC_LIBRARY            : ${MKL_USE_SINGLE_DYNAMIC_LIBRARY}" )
    message("MKL_USE_STATIC_LIBS                       : ${MKL_USE_STATIC_LIBS}" )
    message("BLAS_LIBRARIES                            : ${BLAS_LIBRARIES}" )
    message("LAPACK_LIBRARIES                          : ${LAPACK_LIBRARIES}" )
    message("========================")
    message("")
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
