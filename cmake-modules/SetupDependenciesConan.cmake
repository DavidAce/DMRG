
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
            list(APPEND FOUND_TARGETS mkl::mkl)
        endif()
    else()
        cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
        if(_host_name MATCHES "travis|TRAVIS|Travis|fv-")
            message(STATUS "Setting OpenBLAS dynamic_arch=False")
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=False)
        elseif(_host_name MATCHES "raken")
            message(STATUS
                    "Setting OpenBLAS dynamic_arch=True"
                    "Remember to set environment variable"
                    "   OPENBLAS_CORETYPE=<microarch>"
                    "before launching the executable")
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=True)
        else(_host_name MATCHES "raken")
            message(STATUS "Setting OpenBLAS dynamic_arch=False")
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS openblas:dynamic_arch=False)
        endif()
        find_package(Fortran REQUIRED)
        list(APPEND FOUND_TARGETS gfortran::gfortran)
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

    conan_cmake_run(
            CONANFILE conanfile.txt
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD_TYPE ${CMAKE_BUILD_TYPE}
            BASIC_SETUP CMAKE_TARGETS
            SETTINGS compiler.cppstd=17
            SETTINGS compiler.libcxx=libstdc++11
            PROFILE_AUTO ALL
            ${DMRG_CONAN_OPTIONS}
            BUILD missing
    )

    if(TARGET CONAN_PKG:Eigen3)
        set(eigen_target CONAN_PKG::Eigen3)
    elseif(TARGET CONAN_PKG::eigen)
        set(eigen_target CONAN_PKG::eigen)
    endif()


    if(TARGET ${eigen_target})
        if(TARGET openmp::openmp)
            target_compile_definitions    (${eigen_target} INTERFACE -DEIGEN_USE_THREADS)
        endif()
        if(TB_EIGEN3_BLAS)
            set(EIGEN3_USING_BLAS ON)
            if(TARGET mkl::mkl)
                message(STATUS "Eigen3 will use MKL")
                target_compile_definitions    (${eigen_target} INTERFACE -DEIGEN_USE_MKL_ALL)
                target_compile_definitions    (${eigen_target} INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
                target_link_libraries         (${eigen_target} INTERFACE mkl::mkl)
            elseif(TARGET blas::blas)
                message(STATUS "Eigen3 will use OpenBLAS")
                target_compile_definitions    (${eigen_target} INTERFACE -DEIGEN_USE_BLAS)
                target_compile_definitions    (${eigen_target} INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
                target_link_libraries         (${eigen_target} INTERFACE CONAN_PKG::openblas)
            endif()
        endif()

        cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
        if(_host_name MATCHES "tetralith|triolith")
            # AVX aligns 32 bytes (AVX512 aligns 64 bytes).
            # When running on Tetralith, with march=native, there can be alignment mismatch
            # in ceres which results in a segfault on free memory.
            # Something like "double free or corruption ..."
            #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith

            ### NOTE August 26 2020 ####
            #
            # Ceres started crashing on Tetralith again using -march=native.
            # Tried to solve this issue once and for all.
            # I've tried the following flags during compilation of DMRG++ and ceres-solver:
            #
            #           -DEIGEN_MALLOC_ALREADY_ALIGNED=[none,0,1]
            #           -DEIGEN_MAX_ALIGN_BYTES=[none,16,32]
            #           -march=[none,native]
            #           -std=[none,c++17]
            #
            # Up until now, [0,16,none,none] has worked but now for some reason it stopped now.
            # I noticed the stdc++=17 flag was not being passed on conan builds, so ceres defaulted to -std=c++14 instead.
            # I fixed this in the conanfile.py of the ceres build. The -download-method=fetch method already had this fixed.
            # When no Eigen flags were passed, and ceres-solver finally built with -std=c++17 the issues vanished.
            # In the end what worked was [none,none,native,c++17] in both DMRG++ and ceres-solver.
            # It is important that the same eigen setup is used in all compilation units, and c++17/c++14 seems to
            # make Eigen infer some of the flags differently. In any case, settinc c++17 and no flags for eigen anywhere
            # lets Eigen do its thing in the same way everywhere.

            #                message(STATUS "Applying special Eigen compile definitions for Tetralith: EIGEN_MAX_ALIGN_BYTES=16")
            #                target_compile_definitions(${eigen_target} INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfault?
            #                target_compile_definitions(${eigen_target} INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
        else()
            #                message(STATUS "Applying special Eigen compile definitions for general machines: EIGEN_MAX_ALIGN_BYTES=16")
            #                target_compile_definitions(${eigen_target} INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=1) # May work to fix CERES segfaults!!!
            #                target_compile_definitions(${eigen_target} INTERFACE EIGEN_MAX_ALIGN_BYTES=32)  # May work to fix CERES segfault?
        endif()

    endif()

    # Needs to be linked last
    if(TARGET openmp::openmp)
        list(APPEND FOUND_TARGETS openmp::openmp)
    endif()
endif()
