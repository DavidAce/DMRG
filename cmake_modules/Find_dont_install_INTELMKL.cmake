

# If you want to use Eigen and MKL you need to add #define EIGEN_USE_MKL_ALL before including <Eigen..>
# Adding a -DEIGEN_USE_MKL_ALL here may conflict with arpack++

set(MKL_USE_STATIC_LIBS OFF)
set(MKL_MULTI_THREADED OFF)
set(MKL_USE_SINGLE_DYNAMIC_LIBRARY OFF)
find_package(MKL)
if (MKL_FOUND)
    #    if(MKL_MULTI_THREADED)
    #        list(APPEND MKL_LIBRARIES -lpthread -lm)
    #    endif()
    add_definitions(-DMKL_AVAILABLE)
#    include(cmake/FindGFortran.cmake)     ### For Fortran library
    # To print all variables, use the code below:
    #
#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endforeach()
#    if(MKL_USE_STATIC_LIBS)
#        set(MKL_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
#    else()
#        set(MKL_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
#    endif()
    set(COMPILE_FLAGS "-DMKL_LP64 -m64 -I${MKL_ROOT}/include")
    set(LINK_FLAGS)
    # Link with lmkl_gf_lp64 so that arpack-ng works ( it is compiled with gnu fortran and not intels fortran)
    if(MKL_MULTI_THREADED)
#        set(ENV{OMP_NUM_THREADS} 8)
#        set(ENV{MKL_NUM_THREADS} 8)
        set(LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
    else()
        set(LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    endif()


        #    add_library(mkl::intel              UNKNOWN IMPORTED)
#    add_library(mkl::intel_thread       UNKNOWN IMPORTED)
#    add_library(mkl::core               UNKNOWN IMPORTED)
#    add_library(mkl::cdft_core          UNKNOWN IMPORTED)
#    add_library(mkl::omp                UNKNOWN IMPORTED)
#    add_library(mkl::gf                 UNKNOWN IMPORTED)


#    LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
#
#    LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
#
#    COMPILE_FLAGS "-DMKL_LP64 -m64 -I${MKL_ROOT}/include"
#    set_target_properties(mkl::intel                PROPERTIES  INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR}  IMPORTED_LOCATION   ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64${MKL_SUFFIX})
#    set_target_properties(mkl::intel_thread         PROPERTIES  INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR}  IMPORTED_LOCATION   ${MKL_ROOT}/lib/intel64/libmkl_intel_thread${MKL_SUFFIX})
#    set_target_properties(mkl::core                 PROPERTIES  INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR}  IMPORTED_LOCATION   ${MKL_ROOT}/lib/intel64/libmkl_core${MKL_SUFFIX})
#    set_target_properties(mkl::cdft_core            PROPERTIES  INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR}  IMPORTED_LOCATION   ${MKL_ROOT}/lib/intel64/libmkl_cdft_core${MKL_SUFFIX})
#    set_target_properties(mkl::gf                   PROPERTIES  INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR}  IMPORTED_LOCATION   ${MKL_ROOT}/lib/intel64/libmkl_gf_lp64${MKL_SUFFIX})
    add_library(mkl INTERFACE)
    set_target_properties(mkl PROPERTIES
#            INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
#            INTERFACE_LINK_LIBRARIES "mkl::intel;mkl::intel_thread;mkl::core;mkl::cdft_core;mkl::gf"
            INTERFACE_LINK_LIBRARIES "${MKL_LIBRARIES}"
#            COMPILE_FLAGS "-DMKL_LP64 -m64 -I${MKL_ROOT}/include"
#            LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
            #            LINK_FLAGS "-L${MKL_ROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
            )


    add_library(blas INTERFACE)
    set_target_properties(blas PROPERTIES INTERFACE_LINK_LIBRARIES mkl)

    add_library(lapack INTERFACE)
    set_target_properties(lapack PROPERTIES INTERFACE_LINK_LIBRARIES mkl)

    target_link_libraries(${PROJECT_NAME} mkl)
    target_link_libraries(${PROJECT_NAME} blas)
    target_link_libraries(${PROJECT_NAME} lapack)
    target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INCLUDE_DIR})

    set(BLAS_LIBRARIES   ${MKL_LIBRARIES})
    set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
#    get_target_property(BLAS_LIBRARIES blas IMPORTED_LOCATION)
#    get_target_property(LAPACK_LIBRARIES lapack IMPORTED_LOCATION)

    message("")
    message("======MKL SUMMARY ======")
    message("MKL_LIBRARIES                             : ${MKL_LIBRARIES}" )
    message("MKL_RT_LIBRARY                            : ${MKL_RT_LIBRARY}" )
    message("MKL_INCLUDE_DIR                           : ${MKL_INCLUDE_DIR}" )
    message("MKL_ROOT                                  : ${MKL_ROOT}" )
    message("MKL_USE_SINGLE_DYNAMIC_LIBRARY            : ${MKL_USE_SINGLE_DYNAMIC_LIBRARY}" )
    message("MKL_USE_STATIC_LIBS                       : ${MKL_USE_STATIC_LIBS}" )
    message("MKL_ROOT                                  : ${MKL_ROOT}" )
    message("BLAS_FOUND                                : ${BLAS_FOUND}" )
    message("BLA_VENDOR                                : ${BLA_VENDOR}" )
    message("BLA_STATIC                                : ${BLA_STATIC}" )
    message("BLAS_COMPILER_FLAGS                       : ${BLAS_COMPILER_FLAGS}" )
    message("BLAS_LINKER_FLAGS                         : ${BLAS_LINKER_FLAGS}" )
    message("BLAS_LIBRARIES                            : ${BLAS_LIBRARIES}" )
    message("BLAS95_LIBRARIES                          : ${BLAS95_LIBRARIES}" )
    message("LAPACK_FOUND                              : ${LAPACK_FOUND}" )
    message("LAPACK_INCLUDE_DIR                        : ${LAPACK_INCLUDE_DIR}" )
    message("LAPACK_DEFINITIONS                        : ${LAPACK_DEFINITIONS}" )
    message("LAPACK_LINKER_FLAGS                       : ${LAPACK_LINKER_FLAGS}" )
    message("LAPACK_LIBRARIES_DIR                      : ${LAPACK_LIBRARIES_DIR}" )
    message("LAPACK_LIBRARIES                          : ${LAPACK_LIBRARIES}" )
    message("LAPACK_USE_FILE                           : ${LAPACK_USE_FILE}" )

    message("")
    target_compile_options(${PROJECT_NAME} PUBLIC ${COMPILE_FLAGS})
    target_link_libraries(${PROJECT_NAME}  ${LINK_FLAGS})
endif()

#include <mkl.h>

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #            if("${_variableName}" MATCHES "HDF5")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #            endif()
#endforeach()