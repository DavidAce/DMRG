
if(DMRG_PACKAGE_MANAGER MATCHES "conan")

    ##################################################################
    ### Find the conan command                                     ###
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

    #  Make sure we use our own find modules
    list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
    # Find required packages
    find_package(Threads REQUIRED)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    find_package(gfortran REQUIRED)

    ##############################################################################
    ###  Optional Intel MKL support. Uses OpenBLAS as fall-back                ###
    ##############################################################################
    if(DMRG_ENABLE_MKL)
        find_package(MKL COMPONENTS blas lapack gf gnu_thread lp64 REQUIRED)  # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.

        # Passing BLAS_LIBRARIES as-is will result in cmake-conan injecting -o= between each element
        # Replacing each ";" with spaces will work until arpack-ng tries to link as is.
        # Instead we should use the generator expression trick to let arpack understand that these
        # are multiple libraries
        string(REPLACE ";" "$<SEMICOLON>" MKL_LIBRARIES "${MKL_LIBRARIES}")
        list(APPEND CONAN_OPTIONS
             OPTIONS arpack-ng:blas=All
             OPTIONS arpack-ng:blas_libraries=${MKL_LIBRARIES}
             OPTIONS arpack-ng:lapack_libraries=${MKL_LIBRARIES}
             )
    else()
        if(OPENBLAS_DYNAMIC_ARCH)
            list(APPEND CONAN_OPTIONS OPTIONS openblas:dynamic_arch=True)
        else()
            list(APPEND CONAN_OPTIONS OPTIONS openblas:dynamic_arch=False)
        endif()
    endif()

    if(BUILD_SHARED_LIBS)
        list(APPEND CONAN_OPTIONS OPTIONS "*:shared=True")
    else()
        list(APPEND CONAN_OPTIONS OPTIONS "*:shared=False")
    endif()

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND CONAN_OPTIONS OPTIONS "ceres-solver:use_glog=False")
    endif()

    # Copy the current compiler flags to conan
    string(TOUPPER "${CMAKE_BUILD_TYPE}" CONAN_BUILD_TYPE)
    set(CONAN_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    message(STATUS "CONAN_CXXFLAGS: ${CONAN_CXXFLAGS}")
    message(STATUS "CONAN_CFLAGS  : ${CONAN_CFLAGS}")
    message(STATUS "CONAN_LDFLAGS : ${CONAN_LDFLAGS}")
    message(STATUS "CONAN_BUILD   : ${CONAN_BUILD}")

    conan_add_remote(CONAN_COMMAND ${CONAN_COMMAND} NAME conan-dmrg INDEX 0 URL https://thinkstation.duckdns.org/artifactory/api/conan/conan-dmrg)
    conan_cmake_autodetect(CONAN_AUTODETECT)
    conan_cmake_install(
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD ${CONAN_BUILD} missing outdated cascade
            UPDATE
            GENERATOR CMakeDeps
            SETTINGS ${CONAN_AUTODETECT}
            INSTALL_FOLDER ${CMAKE_BINARY_DIR}/conan
            ENV libunwind:LDFLAGS=-fcommon
            ENV libunwind:CXXFLAGS=-fcommon
            ENV libunwind:CFLAGS=-fcommon
            ENV CC=${CMAKE_C_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
            ENV CXX=${CMAKE_CXX_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
            ENV CXXFLAGS=${CONAN_CXXFLAGS}
            ENV CFLAGS=${CONAN_CFLAGS}
            ENV LDFLAGS=${CONAN_LDFLAGS}
            ENV VERBOSE=1
            ${CONAN_OPTIONS}
            PATH_OR_REFERENCE ${CMAKE_SOURCE_DIR}
    )

    ##################################################################
    ### Find all the things!                                       ###
    ##################################################################
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/conan)
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/conan)
    list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    # Use CONFIG to avoid MODULE mode. This is recommended for the cmake_find_package_multi generator

    find_package(CLI11 2.2.0 REQUIRED CONFIG)
    find_package(Eigen3 3.4 REQUIRED CONFIG)
    find_package(h5pp 1.10.0 REQUIRED CONFIG)
    find_package(fmt 8.1.1 REQUIRED CONFIG)
    find_package(spdlog 1.10.0 REQUIRED CONFIG)
    find_package(arpack-ng 3.8.0 REQUIRED CONFIG)
    find_package(arpack++ 2.3.0 REQUIRED CONFIG)
    find_package(Ceres 2.1.0 REQUIRED CONFIG)
    find_package(Backward 1.6 REQUIRED CONFIG)
    if(NOT TARGET arpack-ng::arpack-ng)
        message(FATAL_ERROR "")
    endif()
    if(NOT DMRG_ENABLE_MKL)
        find_package(OpenBLAS 0.3.20 REQUIRED CONFIG)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE OPENBLAS_AVAILABLE)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_float=std::complex<float>)
        target_compile_definitions(OpenBLAS::OpenBLAS INTERFACE lapack_complex_double=std::complex<double>)
        #For convenience, define these targes
        add_library(BLAS::BLAS ALIAS OpenBLAS::OpenBLAS)
        add_library(LAPACK::LAPACK ALIAS OpenBLAS::OpenBLAS)
        add_library(lapacke::lapacke ALIAS OpenBLAS::OpenBLAS)
    endif()

    target_link_libraries(dmrg-deps INTERFACE
                          CLI11::CLI11
                          h5pp::h5pp
                          arpack++::arpack++
                          primme::primme
                          Ceres::ceres
                          Backward::Backward
                          BLAS::BLAS
                          )

    # Fix issue with Ceres linking to cuda
    find_package(CUDA QUIET) # Same call as when building Ceres
    if(CUDA_FOUND)
        message("-- Found CUDA version ${CUDA_VERSION}: "
                "${CUDA_LIBRARIES};"
                "${CUDA_cusolver_LIBRARY};"
                "${CUDA_cusparse_LIBRARY};"
                "${CUDA_CUBLAS_LIBRARIES}"
                )
        target_link_libraries(dmrg-deps INTERFACE ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
    else()
        target_compile_definitions(dmrg-deps INTERFACE CERES_NO_CUDA)
    endif()

endif()
