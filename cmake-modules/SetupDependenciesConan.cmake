
if(DMRG_DOWNLOAD_METHOD MATCHES "conan")
    #  Make sure we use DMRG's own find modules
    list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake-modules)
    ##############################################################################
    ###  Required OpenMP support                                               ###
    ###  Note that Clang has some  trouble with static openmp and that         ###
    ###  and that static openmp is not recommended. This tries to enable       ###
    ###  static openmp anyway because I find it useful. Installing             ###
    ###  libiomp5 might help for shared linking.                               ###
    ##############################################################################
    find_package(OpenMP REQUIRED) # Uses DMRG's own find module


    ##################################################################
    ### Install conan-modules/conanfile.txt dependencies          ###
    ### This uses conan to get spdlog,eigen3,h5pp,ceres-solver    ###
    ###    ceres-solver/2.0.0@davidace/development                ###
    ###    h5pp/1.7.3@davidace/stable                             ###
    ###    eigen/3.3.7@davidace/patched                           ###
    ##################################################################

    if(DMRG_ENABLE_MKL)
        find_package(Fortran REQUIRED)
        include(cmake-modules/SetupMKL.cmake)         # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
        if(TARGET mkl::mkl)
            include(cmake-modules/getExpandedTarget.cmake)
            expand_target_libs(mkl::mkl MKL_LIBRARIES)
            # Passing MKL_LIBRARIES as-is will result in cmake-conan injecting -o= between each element
            # Replacing each ";" with spaces will work until arpack-ng tries to link as is.
            # Instead we should use the generator expression trick to let arpack understand that these
            # are multiple libraries
            string (REPLACE ";" "$<SEMICOLON>" MKL_LIBRARIES "${MKL_LIBRARIES}")
            list(APPEND DMRG_CONAN_OPTIONS
                    OPTIONS arpack-ng:blas=All
                    OPTIONS arpack-ng:blas_libraries=${MKL_LIBRARIES}
                    OPTIONS arpack-ng:lapack_libraries=${MKL_LIBRARIES}
                    OPTIONS ceres-solver:blas=All
                    OPTIONS ceres-solver:blas_libraries=${MKL_LIBRARIES}
                    OPTIONS ceres-solver:lapack_libraries=${MKL_LIBRARIES}
                    )
            list(APPEND NATIVE_TARGETS mkl::mkl)
        endif()
    else()
        cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
        if(${_host_name} MATCHES "travis|TRAVIS|Travis|fv-")
            message(STATUS "Setting dynamic arch for openblas")
            list(APPEND DMRG_CONAN_OPTIONS
                    OPTIONS openblas:dynamic_arch=True)
        endif()
        find_package(Fortran REQUIRED)
        list(APPEND NATIVE_TARGETS gfortran::gfortran)
    endif()


    find_program (
            CONAN_COMMAND
            conan
            HINTS ${CONAN_PREFIX} $ENV{CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONDA_PREFIX}
            PATHS $ENV{HOME}/anaconda3  $ENV{HOME}/miniconda3 $ENV{HOME}/anaconda $ENV{HOME}/miniconda $ENV{HOME}/.conda
            PATH_SUFFIXES bin envs/dmrg/bin
    )
    if(NOT CONAN_COMMAND)
        message(FATAL_ERROR "Could not find conan program executable")
    else()
        message(STATUS "Found conan: ${CONAN_COMMAND}")
    endif()

    # Download cmake-conan automatically, you can also just copy the conan.cmake file
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()

    include(${CMAKE_BINARY_DIR}/conan.cmake)
    conan_add_remote(NAME conan-center       URL https://conan.bintray.com)
    conan_add_remote(NAME conan-community    URL https://api.bintray.com/conan/conan-community/conan)
    conan_add_remote(NAME bincrafters        URL https://api.bintray.com/conan/bincrafters/public-conan)
    conan_add_remote(NAME conan-dmrg INDEX 1 URL https://api.bintray.com/conan/davidace/conan-dmrg)

    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        # Let it autodetect libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # There is no libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        list(APPEND DMRG_CONAN_SETTINGS SETTINGS compiler.libcxx=libstdc++11)
    endif()
    conan_cmake_run(
            CONANFILE conanfile.txt
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD_TYPE ${CMAKE_BUILD_TYPE}
            BASIC_SETUP CMAKE_TARGETS
            SETTINGS compiler.cppstd=17
            SETTINGS compiler.libcxx=libstdc++11
            ${DMRG_CONAN_SETTINGS}
            ${DMRG_CONAN_OPTIONS}
            BUILD missing
    )



    if(TARGET CONAN_PKG::Eigen3 AND TARGET openmp::openmp)
        target_compile_definitions    (CONAN_PKG::Eigen3 INTERFACE -DEIGEN_USE_THREADS)
    endif()
    if(TARGET CONAN_PKG::Eigen3)
        if(TARGET mkl::mkl)
            message(STATUS "Eigen3 will use MKL")
            set(EIGEN3_USING_BLAS ON)
            target_compile_definitions    (CONAN_PKG::Eigen3 INTERFACE -DEIGEN_USE_MKL_ALL)
            target_compile_definitions    (CONAN_PKG::Eigen3 INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (CONAN_PKG::Eigen3 INTERFACE mkl::mkl)
        else ()
            message(STATUS "Eigen3 will use OpenBLAS")
            set(EIGEN3_USING_BLAS ON)
            target_compile_definitions    (CONAN_PKG::Eigen3 INTERFACE -DEIGEN_USE_BLAS)
            target_compile_definitions    (CONAN_PKG::Eigen3 INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (CONAN_PKG::Eigen3 INTERFACE  CONAN_PKG::openblas)
        endif()


            cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
            if(_host_name MATCHES "tetralith|triolith")
                # AVX aligns 32 bytes (AVX512 aligns 64 bytes).
                # When running on Tetralith, with march=native, there can be alignment mismatch
                # in ceres which results in a segfault on free memory.
                # Something like "double free or corruption ..."
                #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith
                message(STATUS "Applying special Eigen compile definitions for Tetralith: EIGEN_MAX_ALIGN_BYTES=16")
                target_compile_definitions(CONAN_PKG::Eigen3INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfaults!!!
                target_compile_definitions(CONAN_PKG::Eigen3 INTERFACE EIGEN_MAX_ALIGN_BYTES=16)
            else()
                message(STATUS "Applying special Eigen compile definitions for general machines: EIGEN_MAX_ALIGN_BYTES=16")
                target_compile_definitions(CONAN_PKG::Eigen3 INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=1) # May work to fix CERES segfaults!!!
                target_compile_definitions(CONAN_PKG::Eigen3 INTERFACE EIGEN_MAX_ALIGN_BYTES=32)
            endif()

    endif()

    if(TARGET openmp::openmp)
        list(APPEND NATIVE_TARGETS openmp::openmp)
    endif()
endif()
