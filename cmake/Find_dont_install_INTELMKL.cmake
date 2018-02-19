

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# Adding a -DEIGEN_USE_MKL_ALL here may conflict with arpack++

set(MKL_USE_STATIC_LIBS OFF)
set(MKL_MULTI_THREADED OFF)
find_package(MKL)
if (MKL_FOUND)
    if(MKL_MULTI_THREADED)
        list(APPEND MKL_LIBRARIES -lpthread -lm)
    endif()
    add_definitions(-DMKL_AVAILABLE)
    include(cmake/FindGFortran.cmake)     ### For Fortran library
    message("MKL_LIBRARIES: ${MKL_LIBRARIES}" )
    message("MKL_INCLUDE  : ${MKL_INCLUDE_DIR}" )
    message("MKL_ROOT     : ${MKL_ROOT}" )
    message("MKLROOT      : ${MKLROOT}" )
    set(BLAS_FIND_QUIETLY ON)
    set(BLAS_DIR ${MKL_ROOT}/lib/intel64)
    set(BLAS_LIBDIR ${MKL_ROOT}/lib/intel64)
    set(BLAS_VERBOSE ON)
    find_package(BLAS REQUIRED)
    set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
    find_package(LAPACK REQUIRED)

    target_include_directories(${PROJECT_NAME} PRIVATE ${MKL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
    target_compile_options(${PROJECT_NAME} PUBLIC -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-variable)                                   ### Common options
endif()



#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #            if("${_variableName}" MATCHES "HDF5")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #            endif()
#endforeach()