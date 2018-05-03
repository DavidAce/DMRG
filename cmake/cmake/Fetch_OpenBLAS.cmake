
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()

message(STATUS "SEARCHING FOR OpenBLAS IN SYSTEM...")
set(BLA_VENDOR OpenBLAS)
set(BLAS_VERBOSE ON)

find_package(BLAS)
if(BLAS_LIBRARIES)
    message(STATUS "OpenBLAS FOUND IN SYSTEM: ${BLAS_LIBRARIES}")
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    set(BLAS_LOCATION ${BLAS_LIBRARIES})

else()
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")

    enable_language(Fortran)
    include(ExternalProject)
    ExternalProject_Add(project_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.2.20
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""

            CMAKE_ARGS
            -j8
            -DBUILD_SHARED_LIBS=OFF
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w
            -DCMAKE_CXX_FLAGS=-w
            -DCMAKE_FORTRAN_FLAGS=-w
            -DCMAKE_BUILD_TYPE=Release
            )

    ExternalProject_Get_Property(project_OpenBLAS INSTALL_DIR)

    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(BLAS_LOCATION ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas project_OpenBLAS)
    add_dependencies(lapack project_OpenBLAS)

endif()


set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION ${BLAS_LOCATION}
        INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)
set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION ${BLAS_LOCATION}
        INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)
target_link_libraries(${PROJECT_NAME} blas)
target_link_libraries(${PROJECT_NAME} lapack)
target_include_directories(${PROJECT_NAME} PRIVATE ${BLAS_INCLUDE_DIRS})
#For convenience, define these variables
get_target_property(BLAS_LIBRARIES blas IMPORTED_LOCATION)
get_target_property(LAPACK_LIBRARIES lapack IMPORTED_LOCATION)
