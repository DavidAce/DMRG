
# If the INTEL MKL library has already been loaded, skip the rest.
if(TARGET mkl::mkl)
    return()
endif()


include(cmake-modules/FindOpenBLAS.cmake)
find_OpenBLAS()

# To print all variables, use the code below:

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     if("${_variableName}" MATCHES "BLAS" OR "${_variableName}" MATCHES "blas" OR "${_variableName}" MATCHES "Blas")
#         message(STATUS "${_variableName}=${${_variableName}}")
#     endif()
# endforeach()

if(TARGET openblas::openblas)
#    message(STATUS "OpenBLAS found")

elseif(NOT ${DMRG_DOWNLOAD_METHOD} MATCHES "none")
    message(STATUS "OpenBLAS will be installed into ${CMAKE_INSTALL_PREFIX}/OpenBLAS")
    set(OpenBLAS_MULTITHREADED 1 )
    if(TARGET openmp::openmp)
        set(OpenBLAS_ENABLE_OPENMP 1) # Openmp doesnt work on clang it seems
    else()
        set(OpenBLAS_ENABLE_OPENMP 0) # Openmp doesnt work on clang it seems
    endif()
    include(cmake-modules/getExpandedTarget.cmake)
    expand_target_libs(gfortran::gfortran GFORTRAN_LIBS)
    list(FILTER GFORTRAN_LIBS INCLUDE REGEX "gfortran|quadmath")
    list(JOIN GFORTRAN_LIBS " " GFORTRAN_LIBS_SPACE_SEP_STRING)
    list(GET  GFORTRAN_LIBS 0 GFORTRAN_LIB)
    get_filename_component(GFORTRAN_PATH ${GFORTRAN_LIB} DIRECTORY)

    set(LDFLAGS "-L${GFORTRAN_PATH} ${GFORTRAN_LIBS_SPACE_SEP_STRING}")
    set(FFLAGS  "-O3 -Wno-maybe-uninitialized -Wno-conversion -Wno-unused-but-set-variable -Wno-unused-variable")
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DGFORTRAN_PATH:PATH=${GFORTRAN_PATH})
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DTARGET:STRING=${OPENBLAS_MARCH})
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DUSE_THREAD:BOOL=${OpenBLAS_MULTITHREADED})
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DUSE_OPENMP:BOOL=${OpenBLAS_ENABLE_OPENMP})
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DLDFLAGS:STRING=${LDFLAGS})
    list(APPEND OpenBLAS_CMAKE_OPTIONS  -DFFLAGS:STRING=${FFLAGS})
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(OpenBLAS "${CMAKE_INSTALL_PREFIX}" "${OpenBLAS_CMAKE_OPTIONS}")
    find_OpenBLAS()
    if(TARGET openblas::openblas)
        message(STATUS "OpenBLAS installed successfully")
    else()
        message(FATAL_ERROR "OpenBLAS could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency OpenBLAS not found and DMRG_DOWNLOAD_METHOD = ${DMRG_DOWNLOAD_METHOD}")
endif()
