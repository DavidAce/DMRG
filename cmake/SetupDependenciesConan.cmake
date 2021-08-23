
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
    endif()

    unset(CONAN_BUILD_INFO)
    unset(CONAN_BUILD_INFO CACHE)
    find_file(CONAN_BUILD_INFO
            conanbuildinfo.cmake
            HINTS ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR}
            NO_DEFAULT_PATH)

    if(CONAN_BUILD_INFO)
        ##################################################################
        ### Use pre-existing conanbuildinfo.cmake                      ###
        ### This avoids having to run conan again                      ###
        ##################################################################
        message(STATUS "Detected Conan build info: ${CONAN_BUILD_INFO}")
        include(${CONAN_BUILD_INFO})
        conan_basic_setup(KEEP_RPATHS TARGETS)
    else()

        ##################################################################
        ### Use cmake-conan integration to launch conan                ###
        ### Install dependencies from conanfile.txt                    ###
        ### This uses conan to get spdlog,eigen3,h5pp,ceres-solver     ###
        ###    ceres-solver/2.0.0@davidace/development                 ###
        ###    h5pp/1.8.5@davidace/stable                              ###
        ###    eigen/3.3.9@davidace/patched                            ###
        ##################################################################
        unset(CONAN_COMMAND CACHE)
        find_program (
                CONAN_COMMAND
                conan
                HINTS ${CONAN_PREFIX} $ENV{CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONDA_PREFIX}
                PATHS
                $ENV{HOME}/anaconda3
                $ENV{HOME}/miniconda3
                $ENV{HOME}/anaconda
                $ENV{HOME}/miniconda
                $ENV{HOME}/.local
                $ENV{HOME}/.conda
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
            file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.16.1/conan.cmake"
                    "${CMAKE_BINARY_DIR}/conan.cmake")
        endif()

        if(BUILD_SHARED_LIBS)
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS "*:shared=True")
        else()
            list(APPEND DMRG_CONAN_OPTIONS OPTIONS "*:shared=False")
        endif()


        include(${CMAKE_BINARY_DIR}/conan.cmake)
        conan_add_remote(NAME conan-dmrg URL http://thinkstation.duckdns.org:8081/artifactory/api/conan/conan-dmrg)

        conan_cmake_run(
                CONANFILE conanfile.txt
                CONAN_COMMAND ${CONAN_COMMAND}
                BUILD_TYPE ${CMAKE_BUILD_TYPE}
                BASIC_SETUP CMAKE_TARGETS
                SETTINGS compiler.cppstd=17
                SETTINGS compiler.libcxx=libstdc++11
                ENV libunwind:LDFLAGS=-fcommon
                ENV libunwind:CXXFLAGS=-fcommon
                ENV libunwind:CFLAGS=-fcommon
                PROFILE_AUTO ALL
                ${DMRG_CONAN_OPTIONS}
                KEEP_RPATHS
                BUILD missing
                BUILD openblas # This builds openblas everytime on github actions
        )

    endif()

    if(TARGET CONAN_PKG::eigen AND DMRG_ENABLE_THREADS)
        target_compile_definitions(CONAN_PKG::eigen INTERFACE EIGEN_USE_THREADS)
    endif()

    if(TARGET CONAN_PKG::eigen)
        if(TARGET mkl::mkl)
            message(STATUS "Eigen3 will use MKL")
            target_compile_definitions    (CONAN_PKG::eigen INTERFACE EIGEN_USE_MKL_ALL)
            target_compile_definitions    (CONAN_PKG::eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (CONAN_PKG::eigen INTERFACE mkl::mkl)
        elseif(TARGET CONAN_PKG::openblas)
            message(STATUS "Eigen3 will use OpenBLAS")
            target_compile_definitions    (CONAN_PKG::eigen INTERFACE EIGEN_USE_BLAS)
            target_compile_definitions    (CONAN_PKG::eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
            target_link_libraries         (CONAN_PKG::eigen INTERFACE CONAN_PKG::openblas)
        endif()

        cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
        if(_host_name MATCHES "tetralith|triolith")
            # AVX aligns 32 bytes (AVX512 aligns 64 bytes).
            # When running on Tetralith, with march=native, there can be alignment mismatch
            # in ceres which results in a segfault on free memory.
            # Something like "double free or corruption ..."
            #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith

            ### NOTE October 4 2020 ####
            #
            # Another flag that seems to fix weird release-only bugs is
            #       -fno-strict-aliasing

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
            # I fixed this in the conanfile.py of the ceres build. The -package-manager=cmake method already had this fixed.
            # When no Eigen flags were passed, and ceres-solver finally built with -std=c++17 the issues vanished.
            # In the end what worked was [none,none,native,c++17] in both DMRG++ and ceres-solver.
            # It is important that the same eigen setup is used in all compilation units, and c++17/c++14 seems to
            # make Eigen infer some of the flags differently. In any case, settinc c++17 and no flags for eigen anywhere
            # lets Eigen do its thing in the same way everywhere.

#            message(STATUS "Applying special Eigen compile definitions for Tetralith: EIGEN_MAX_ALIGN_BYTES=16")
#            target_compile_definitions(CONAN_PKG::eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfault?
#            target_compile_definitions(CONAN_PKG::eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
        else()
#            message(STATUS "Applying special Eigen compile definitions for general machines: EIGEN_MAX_ALIGN_BYTES=16")
#            target_compile_definitions(CONAN_PKG::eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfaults!!!
#            target_compile_definitions(CONAN_PKG::eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
        endif()
    else()
        message(FATAL_ERROR "Target not defined: CONAN_PKG::eigen")
    endif()


    # Make aliases
    add_library(cxxopts::cxxopts    ALIAS CONAN_PKG::cxxopts)
    add_library(Eigen3::Eigen       ALIAS CONAN_PKG::eigen)
    add_library(h5pp::h5pp          ALIAS CONAN_PKG::h5pp)
    add_library(fmt::fmt            ALIAS CONAN_PKG::fmt)
    add_library(spdlog::spdlog      ALIAS CONAN_PKG::spdlog)
    add_library(arpack::arpack++    ALIAS CONAN_PKG::arpack++)
    add_library(Ceres::ceres        ALIAS CONAN_PKG::ceres-solver)
    add_library(unwind::unwind      ALIAS CONAN_PKG::libunwind)
    if(TARGET CONAN_PKG::openblas)
        add_library(OpenBLAS::OpenBLAS  ALIAS CONAN_PKG::openblas)
        #For convenience, define these targes
        add_library(BLAS::BLAS ALIAS CONAN_PKG::openblas)
        add_library(LAPACK::LAPACK ALIAS CONAN_PKG::openblas)
        add_library(lapacke::lapacke  ALIAS CONAN_PKG::openblas)

    endif()
endif()
