

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# Adding a -DEIGEN_USE_MKL_ALL here may conflict with arpack++

set(MKL_USE_STATIC_LIBS ON)
set(MKL_MULTI_THREADED OFF)
set(MKL_USE_SINGLE_DYNAMIC_LIBRARY OFF)
find_package(MKL)
if (MKL_FOUND)
    add_definitions(-DMKL_AVAILABLE)
    set(MKL_FLAGS -DMKL_LP64 -m64 -I${MKL_ROOT}/include)
    # Link with lmkl_gf_lp64 so that arpack-ng works (it is compiled with gnu fortran and not intels fortran)
    if(MKL_MULTI_THREADED)
        set(MKL_LFLAGS -L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl)
    else()
        set(MKL_LFLAGS -L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl)
    endif()

    add_library(mkl INTERFACE)
    set_target_properties(mkl PROPERTIES
            INTERFACE_LINK_LIBRARIES "${MKL_LIBRARIES}"
            INTERFACE_COMPILE_OPTIONS "${MKL_FLAGS}"
            INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
            INTERFACE_POSITION_INDEPENDENT_CODE ON
            )
    target_link_libraries(mkl INTERFACE ${MKL_LFLAGS})

    target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PUBLIC mkl)
    target_compile_options(${PROJECT_NAME} PUBLIC ${MKL_FLAGS})



    # Make the rest of the build structure aware of blas and lapack included in MKL.

    add_library(blas INTERFACE)
    set_target_properties(blas PROPERTIES
            INTERFACE_LINK_LIBRARIES mkl)

    add_library(lapack INTERFACE)
    set_target_properties(lapack PROPERTIES
            INTERFACE_LINK_LIBRARIES mkl)


    set(BLAS_LIBRARIES   ${MKL_LIBRARIES})
    set(LAPACK_LIBRARIES ${MKL_LIBRARIES})



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
    message("MKL_USE_SINGLE_DYNAMIC_LIBRARY            : ${MKL_USE_SINGLE_DYNAMIC_LIBRARY}" )
    message("MKL_USE_STATIC_LIBS                       : ${MKL_USE_STATIC_LIBS}" )
    message("BLAS_LIBRARIES                            : ${BLAS_LIBRARIES}" )
    message("LAPACK_LIBRARIES                          : ${LAPACK_LIBRARIES}" )
    message("========================")
    message("")
#    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES} ${MKL_LFLAGS})
endif()


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #            if("${_variableName}" MATCHES "HDF5")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #            endif()
#endforeach()