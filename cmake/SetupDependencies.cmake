
# Find threads
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Threads::Threads INTERFACE rt dl)

# Find OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)
target_link_libraries(flags INTERFACE OpenMP::OpenMP_CXX)


# Setup dependencies
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)

# Install dependencies that are not in conan.
include(cmake/InstallPackage.cmake)
install_package(primme MODULE)

# Link dependencies
if(NOT TARGET deps)
    add_library(deps INTERFACE)
endif()
target_link_libraries(deps INTERFACE
        CLI11::CLI11
        h5pp::h5pp
        arpack++::arpack++
        primme::primme
        Ceres::ceres
        BLAS::BLAS
        )
if(TARGET unwind::unwind)
    target_link_libraries(deps INTERFACE unwind::unwind)
    target_compile_definitions(deps INTERFACE DMRG_HAS_UNWIND=1)
endif()

if(TARGET Backward::Backward)
    target_link_libraries(deps INTERFACE Backward::Backward)
endif()


# Configure Eigen
if(TARGET Eigen3::Eigen)
    target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    if(TARGET mkl::mkl)
        message(STATUS "Eigen3 will use MKL")
        target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
        target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
        target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
    elseif(TARGET OpenBLAS::OpenBLAS)
        message(STATUS "Eigen3 will use OpenBLAS")
        target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
        target_compile_definitions    (Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
        target_link_libraries         (Eigen3::Eigen INTERFACE OpenBLAS::OpenBLAS)
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
        #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfault?
        #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
    else()
        #            message(STATUS "Applying special Eigen compile definitions for general machines: EIGEN_MAX_ALIGN_BYTES=16")
        #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfaults!!!
        #            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)  # May work to fix CERES segfault?
    endif()
else()
    message(FATAL_ERROR "Target not defined: Eigen3::Eigen")
endif()

