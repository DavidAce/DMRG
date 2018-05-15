
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

message(STATUS "SEARCHING FOR BLAS IN SYSTEM...")
set(BLA_VENDOR Open)
set(BLAS_VERBOSE ON)

find_package(BLAS)

if(BLAS_FOUND)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    message(STATUS "SEARCHING FOR LAPACK IN SYSTEM...")
    set(BLAS_DIR ${BLAS_openblas_LIBRARY}) # Let Lapack find the same BLAS implementation
    set(BLA_VENDOR OpenBLAS)        # BLA_VENDOR variable is different for this cmake module
    find_package(LAPACK)
    if(LAPACK_FOUND)
        message(STATUS "LAPACK FOUND IN SYSTEM: ${LAPACK_openblas_LIBRARY}")
        add_library(blas UNKNOWN IMPORTED)
        add_library(lapack UNKNOWN IMPORTED)
        set(BLAS_LOCATION ${BLAS_openblas_LIBRARY})
        set(LAPACK_LOCATION ${LAPACK_openblas_LIBRARY})
    endif()
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        if("${_variableName}" MATCHES "blas"  OR "${_variableName}" MATCHES "BLAS" OR "${_variableName}" MATCHES "LAPACK")
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()
endif()
#exit (1)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")

    enable_language(Fortran)
    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.2.20
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
#            UPDATE_COMMAND ""
#            TEST_COMMAND ""
#
#            CMAKE_ARGS
#            -j8
#            -DBUILD_SHARED_LIBS=OFF
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
#            -DCMAKE_C_FLAGS=-w
#            -DCMAKE_CXX_FLAGS=-w
#            -DCMAKE_FORTRAN_FLAGS=-w
#            -DCMAKE_BUILD_TYPE=Release



            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) USE_THREAD=0 USE_OPENMP=0 NO_LAPACKE=1 NO_CBLAS=1 BINARY64=1
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            )

    ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)

    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(BLAS_LOCATION ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LOCATION ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas library_OpenBLAS)
    add_dependencies(lapack library_OpenBLAS)

endif()


set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION ${BLAS_LOCATION}
        INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)
set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION ${LAPACK_LOCATION}
        INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)
target_link_libraries(${PROJECT_NAME} blas -lpthread)
target_link_libraries(${PROJECT_NAME} lapack)
target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})
#For convenience, define these variables
get_target_property(BLAS_LIBRARIES blas IMPORTED_LOCATION)
get_target_property(LAPACK_LIBRARIES lapack IMPORTED_LOCATION)
