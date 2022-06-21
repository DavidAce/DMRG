
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
        expand_target_libs(mkl::mkl MKL_LIBRARIES)
        # Passing BLAS_LIBRARIES as-is will result in cmake-conan injecting -o= between each element
        # Replacing each ";" with spaces will work until arpack-ng tries to link as is.
        # Instead we should use the generator expression trick to let arpack understand that these
        # are multiple libraries
        string(REPLACE ";" "$<SEMICOLON>" MKL_LIBRARIES "${MKL_LIBRARIES}")
        list(APPEND DMRG_CONAN_OPTIONS
                OPTIONS arpack-ng:blas=All
                OPTIONS arpack-ng:blas_libraries=${MKL_LIBRARIES}
                OPTIONS arpack-ng:lapack_libraries=${MKL_LIBRARIES}
                )
    else()
        if(OPENBLAS_DYNAMIC_ARCH)
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=True)
        else()
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=False)
        endif()
    endif()


    ##################################################################
    ### Install dependencies from conanfile.txt                    ###
    ##################################################################
    unset(CONAN_COMMAND CACHE)
    find_program(CONAN_COMMAND conan
            HINTS ${DMRG_CONAN_HINTS}
            PATH_SUFFIXES ${DMRG_CONAN_PATH_SUFFIXES}
            REQUIRED)


    # Download cmake-conan integrator
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan/conan.cmake"
                EXPECTED_HASH MD5=81d5eab13a49f43527e35a90bfac6960
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

    conan_add_remote(CONAN_COMMAND ${CONAN_COMMAND} NAME conan-dmrg INDEX 0 URL https://thinkstation.duckdns.org/artifactory/api/conan/conan-dmrg)
    conan_cmake_autodetect(CONAN_AUTODETECT)
    conan_cmake_install(
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD missing outdated cascade
            GENERATOR CMakeDeps
            SETTINGS ${CONAN_AUTODETECT}
            INSTALL_FOLDER ${CMAKE_BINARY_DIR}/conan
            ENV libunwind:LDFLAGS=-fcommon
            ENV libunwind:CXXFLAGS=-fcommon
            ENV libunwind:CFLAGS=-fcommon
            ENV CC=${CMAKE_C_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
            ENV CXX=${CMAKE_CXX_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
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

    find_package(CLI11 2.2.0 REQUIRED CONFIG)
    find_package(Eigen 3.4 REQUIRED CONFIG)
    find_package(h5pp 1.10.0 REQUIRED CONFIG)
    find_package(fmt 8.1.1 REQUIRED CONFIG)
    find_package(spdlog 1.10.0 REQUIRED CONFIG)
    find_package(arpack++ 2.3.0 REQUIRED CONFIG)
    find_package(ceres-solver 2.0.0 REQUIRED CONFIG)
    find_package(libunwind 1.6.2 REQUIRED CONFIG)
    find_package(backward-cpp 1.6 REQUIRED CONFIG)
    if (NOT DMRG_ENABLE_MKL)
        find_package(OpenBLAS 0.3.17 REQUIRED CONFIG)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
        #For convenience, define these targes
        add_library(BLAS::BLAS ALIAS OpenBLAS::OpenBLAS)
        add_library(LAPACK::LAPACK ALIAS OpenBLAS::OpenBLAS)
        add_library(lapacke::lapacke ALIAS OpenBLAS::OpenBLAS)
    endif ()
    if (TARGET eigen::eigen AND NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen ALIAS eigen::eigen)
    endif ()
    if (TARGET ceres-solver::ceres-solver AND NOT TARGET Ceres::Ceres)
        add_library(Ceres::ceres ALIAS ceres-solver::ceres-solver)
    endif ()
    if (TARGET backward-cpp::backward-cpp AND NOT TARGET Backward::Backward)
        add_library(Backward::Backward ALIAS backward-cpp::backward-cpp)
    endif ()
    target_link_libraries(dmrg-deps INTERFACE
            CLI11::CLI11
            h5pp::h5pp
            arpack++::arpack++
            primme::primme
            ceres-solver::ceres-solver
            BLAS::BLAS
            backward-cpp::backward-cpp
            )

    if (TARGET libunwind::libunwind)
        target_compile_definitions(dmrg-deps INTERFACE DMRG_HAS_UNWIND=1)
        target_link_libraries(dmrg-deps INTERFACE libunwind::libunwind)
    endif ()

    # Configure Eigen
    if (TARGET eigen::eigen)
        target_compile_definitions(eigen::eigen INTERFACE EIGEN_USE_THREADS)
        if (TARGET BLAS::BLAS)
            target_link_libraries(eigen::eigen INTERFACE BLAS::BLAS)
            target_compile_definitions(eigen::eigen INTERFACE EIGEN_USE_BLAS)
            target_compile_definitions(eigen::eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            if (TARGET mkl::mkl)
                message(STATUS "Eigen3 will use MKL")
                target_compile_definitions(eigen::eigen INTERFACE EIGEN_USE_MKL_ALL)
            else ()
                message(STATUS "Eigen3 will use OpenBLAS")
            endif ()
        endif ()
    else ()
        message(FATAL_ERROR "Target not defined: eigen::eigen")
    endif ()
endif()
