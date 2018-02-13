


set(MKL_USE_STATIC_LIBS ON)
set(MKL_MULTI_THREADED ON)
find_package(MKL)
if (MKL_FOUND)
    if(MKL_MULTI_THREADED)
        list(APPEND MKL_LIBRARIES -lpthread -lm)
    endif()


    message("MKL_LIBRARIES: ${MKL_LIBRARIES}" )
    message("MKL_INCLUDE  : ${MKL_INCLUDE_DIR}" )
    message("MKL_ROOT     : ${MKL_ROOT}" )
    message("MKLROOT      : ${MKLROOT}" )
    set(BLAS_DIR ${MKL_ROOT}/lib/intel64)
    set(BLAS_LIBDIR ${MKL_ROOT}/lib/intel64)
    set(BLAS_VERBOSE ON)
    find_package(BLAS REQUIRED)
    set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
    find_package(LAPACK REQUIRED)

    target_include_directories(${PROJECT_NAME} PRIVATE ${MKL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
    target_compile_options(${PROJECT_NAME} PUBLIC -DEIGEN_USE_MKL_ALL -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-variable)                                   ### Common options
endif()



#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    #            if("${_variableName}" MATCHES "HDF5")
#    message(STATUS "${_variableName}=${${_variableName}}")
#    #            endif()
#endforeach()