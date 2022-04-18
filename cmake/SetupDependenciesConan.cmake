
if(DMRG_PACKAGE_MANAGER MATCHES "conan")

    #  Make sure we use DMRG's own find modules
    list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake)
    find_package(Threads REQUIRED)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    find_package(Fortran REQUIRED)
    # Find packages or install if missing

    ##############################################################################
    ###  Optional Intel MKL support. Uses OpenBLAS as fall-back                ###
    ##############################################################################
    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS blas lapack gf gnu_thread lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
        add_library(lapacke::lapacke ALIAS mkl::mkl) # Lapacke is included

        include(cmake/getExpandedTarget.cmake)
        expand_target_libs(BLAS::BLAS BLAS_LIBRARIES)
        # Passing BLAS_LIBRARIES as-is will result in cmake-conan injecting -o= between each element
        # Replacing each ";" with spaces will work until arpack-ng tries to link as is.
        # Instead we should use the generator expression trick to let arpack understand that these
        # are multiple libraries
        string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES "${BLAS_LIBRARIES}")
        list(APPEND DMRG_CONAN_OPTIONS
                OPTIONS arpack-ng:blas=All
                OPTIONS arpack-ng:blas_libraries=${BLAS_LIBRARIES}
                OPTIONS arpack-ng:lapack_libraries=${BLAS_LIBRARIES}
                )
    else()
        if(NOT OPENBLAS_DYNAMIC_ARCH)
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=False)
        endif()
    endif()


    ##################################################################
    ### Install dependencies from conanfile.txt                    ###
    ##################################################################
    unset(CONAN_COMMAND CACHE)
    find_program (CONAN_COMMAND conan
            HINTS ${DMRG_CONAN_HINTS}
            PATH_SUFFIXES ${DMRG_CONAN_PATH_SUFFIXES})
    if(NOT CONAN_COMMAND)
        message(FATAL_ERROR "Could not find conan program executable")
    else()
        message(STATUS "Found conan: ${CONAN_COMMAND}")
    endif()

    # Download cmake-conan integrator
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/release/0.17/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan/conan.cmake"
                EXPECTED_HASH MD5=52a255a933397fdce3d0937f9c737e98
                TLS_VERIFY ON)
    endif()
    include(${CMAKE_BINARY_DIR}/conan/conan.cmake)

    if(BUILD_SHARED_LIBS)
        list(APPEND DMRG_CONAN_OPTIONS OPTIONS "*:shared=True")
    else()
        list(APPEND DMRG_CONAN_OPTIONS OPTIONS "*:shared=False")
    endif()

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND DMRG_CONAN_OPTIONS OPTIONS "ceres-solver:use_glog=False")
    endif()

    conan_add_remote(CONAN_COMMAND ${CONAN_COMMAND} NAME conan-dmrg URL https://thinkstation.duckdns.org/artifactory/api/conan/conan-dmrg)
    conan_cmake_autodetect(CONAN_AUTODETECT)
    conan_cmake_install(
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD missing outdated cascade
            GENERATOR cmake_find_package_multi
            SETTINGS ${CONAN_AUTODETECT}
            INSTALL_FOLDER ${CMAKE_BINARY_DIR}/conan
            ENV libunwind:LDFLAGS=-fcommon
            ENV libunwind:CXXFLAGS=-fcommon
            ENV libunwind:CFLAGS=-fcommon
            ${DMRG_CONAN_OPTIONS}
            PATH_OR_REFERENCE ${CMAKE_SOURCE_DIR}
    )

    ##################################################################
    ### Find all the things!                                       ###
    ##################################################################
    if(NOT CONAN_CMAKE_SILENT_OUTPUT)
        set(CONAN_CMAKE_SILENT_OUTPUT OFF) # Default is off
    endif()
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/conan)
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/conan)
    # Use CONFIG to avoid MODULE mode. This is recommended for the cmake_find_package_multi generator

    find_package(CLI11        2.2.0  REQUIRED CONFIG)
    find_package(Eigen3       3.4    REQUIRED CONFIG)
    find_package(h5pp         1.10.0 REQUIRED CONFIG)
    find_package(fmt          8.1.1  REQUIRED CONFIG)
    find_package(spdlog       1.10.0  REQUIRED CONFIG)
    find_package(arpack++     2.3.0  REQUIRED CONFIG)
    find_package(Ceres        2.0.0  REQUIRED CONFIG)
    find_package(libunwind    1.6.2  REQUIRED CONFIG)
    find_package(Backward     1.6    REQUIRED CONFIG)
    if(NOT DMRG_ENABLE_MKL)
        find_package(OpenBLAS 0.3.17 REQUIRED CONFIG)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        #For convenience, define these targes
        add_library(BLAS::BLAS ALIAS OpenBLAS::OpenBLAS)
        add_library(LAPACK::LAPACK ALIAS OpenBLAS::OpenBLAS)
        add_library(lapacke::lapacke  ALIAS OpenBLAS::OpenBLAS)
    endif()

endif()
