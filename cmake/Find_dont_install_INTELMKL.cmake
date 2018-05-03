

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# Adding a -DEIGEN_USE_MKL_ALL here may conflict with arpack++

set(MKL_USE_STATIC_LIBS OFF)
set(MKL_MULTI_THREADED OFF)
set(MKL_USE_SINGLE_DYNAMIC_LIBRARY ON)
find_package(MKL)
if (MKL_FOUND)
    #    if(MKL_MULTI_THREADED)
    #        list(APPEND MKL_LIBRARIES -lpthread -lm)
    #    endif()
    add_definitions(-DMKL_AVAILABLE)
    include(cmake/FindGFortran.cmake)     ### For Fortran library
    # To print all variables, use the code below:
    #
#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endforeach()

    add_library(mkl SHARED IMPORTED)
#    add_library(blas SHARED IMPORTED)
#    add_library(lapack SHARED IMPORTED)
    set_target_properties(mkl PROPERTIES
            IMPORTED_LOCATION ${MKL_RT_LIBRARY}
            INCLUDE_DIRECTORIES MKL_INCLUDE_DIR)
#    set_target_properties(blas PROPERTIES
#            IMPORTED_LOCATION ${MKL_RT_LIBRARY}
#            INCLUDE_DIRECTORIES MKL_INCLUDE_DIR)
#    set_target_properties(lapack PROPERTIES
#            IMPORTED_LOCATION ${MKL_RT_LIBRARY}
#            INCLUDE_DIRECTORIES MKL_INCLUDE_DIR)

    target_link_libraries(${PROJECT_NAME} mkl)
#    target_link_libraries(${PROJECT_NAME} blas)
#    target_link_libraries(${PROJECT_NAME} lapack)
    target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INCLUDE_DIR})

#    get_target_property(BLAS_LIBRARIES blas IMPORTED_LOCATION)
#    get_target_property(LAPACK_LIBRARIES lapack IMPORTED_LOCATION)

    message("")
    message("======MKL SUMMARY ======")
    message("MKL_LIBRARY         : ${MKL_RT_LIBRARY}" )
    message("MKL_INCLUDE         : ${MKL_INCLUDE_DIR}" )
    message("MKL_ROOT            : ${MKL_ROOT}" )
    message("BLAS_FOUND          : ${BLAS_FOUND}" )
    message("BLA_VENDOR          : ${BLA_VENDOR}" )
    message("BLA_STATIC          : ${BLA_STATIC}" )
    message("BLAS_COMPILER_FLAGS : ${BLAS_COMPILER_FLAGS}" )
    message("BLAS_LINKER_FLAGS   : ${BLAS_LINKER_FLAGS}" )
    message("BLAS_LIBRARIES      : ${BLAS_LIBRARIES}" )
    message("BLAS95_LIBRARIES    : ${BLAS95_LIBRARIES}" )
    message("LAPACK_FOUND        : ${LAPACK_FOUND}" )
    message("LAPACK_INCLUDE_DIR  : ${LAPACK_INCLUDE_DIR}" )
    message("LAPACK_DEFINITIONS  : ${LAPACK_DEFINITIONS}" )
    message("LAPACK_LINKER_FLAGS : ${LAPACK_LINKER_FLAGS}" )
    message("LAPACK_LIBRARIES_DIR: ${LAPACK_LIBRARIES_DIR}" )
    message("LAPACK_LIBRARIES    : ${LAPACK_LIBRARIES}" )
    message("LAPACK_USE_FILE     : ${LAPACK_USE_FILE}" )

    message("")
    target_compile_options(${PROJECT_NAME} PUBLIC  -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-variable)                                   ### Common options
    target_link_libraries(${PROJECT_NAME} -L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl)
endif()

#include <mkl.h>

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #            if("${_variableName}" MATCHES "HDF5")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #            endif()
#endforeach()