
# If the INTEL MKL library has already been loaded, skip the rest.
if(MKL_FOUND)
    return()
endif()


message(STATUS "SEARCHING FOR BLAS IN SYSTEM...")
set(BLA_VENDOR Open)
set(BLAS_VERBOSE OFF)
find_package(BLAS)
if(BLAS_FOUND)
    message(STATUS "BLAS FOUND IN SYSTEM: ${BLAS_openblas_LIBRARY}")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-lopenblas _support_lopenblas)
    check_cxx_compiler_flag(-llapack _support_llapack)
    if(_support_lopenblas AND _support_llapack)


        target_link_libraries(${PROJECT_NAME} openblas)
        target_link_libraries(${PROJECT_NAME} lapack)
        # Make dummy library blas and lapack pointing to openblas
        add_library(blas INTERFACE)
        add_library(lapack INTERFACE)
        set_target_properties(blas PROPERTIES
                INTERFACE_LINK_LIBRARIES "openblas")
        set_target_properties(lapack PROPERTIES
                INTERFACE_LINK_LIBRARIES "lapack;openblas")

        #For convenience, define these variables
        set(BLAS_LIBRARIES     ${BLAS_openblas_LIBRARY})
        set(LAPACK_LIBRARIES   ${BLAS_openblas_LIBRARY})
        return()
    else()
        unset(BLAS_FOUND)
    endif()
endif()
#exit (1)

if(NOT BLAS_FOUND)
    message(STATUS "OpenBLAS will be installed into ${INSTALL_DIRECTORY}/OpenBLAS on first build.")

    enable_language(Fortran)
    include(ExternalProject)
    ExternalProject_Add(library_OpenBLAS
            GIT_REPOSITORY      https://github.com/xianyi/OpenBLAS.git
            GIT_TAG             v0.2.20
            PREFIX              "${INSTALL_DIRECTORY}/OpenBLAS"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND $(MAKE) USE_THREAD=0 USE_OPENMP=0 NO_LAPACKE=1 NO_CBLAS=1 BINARY64=1
            INSTALL_COMMAND $(MAKE) PREFIX=<INSTALL_DIR> install
            )

    ExternalProject_Get_Property(library_OpenBLAS INSTALL_DIR)
    set(BLAS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(BLAS_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(LAPACK_LIBRARIES ${INSTALL_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    add_library(blas UNKNOWN IMPORTED)
    add_library(lapack UNKNOWN IMPORTED)
    add_dependencies(blas library_OpenBLAS)
    add_dependencies(lapack library_OpenBLAS)
    set_target_properties(blas PROPERTIES
            IMPORTED_LOCATION ${BLAS_LIBRARIES}
            INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)
    set_target_properties(lapack PROPERTIES
            IMPORTED_LOCATION ${LAPACK_LIBRARIES}
            INCLUDE_DIRECTORIES BLAS_INCLUDE_DIRS)

    target_link_libraries(${PROJECT_NAME} blas)
    target_link_libraries(${PROJECT_NAME} lapack)
    target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})
endif()



