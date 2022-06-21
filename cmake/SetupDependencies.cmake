
# Find threads
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Threads::Threads INTERFACE rt dl)

# Find OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)
target_link_libraries(dmrg-flags INTERFACE OpenMP::OpenMP_CXX)


# Link all dependencies to dmrg-deps
if (NOT TARGET dmrg-deps)
    add_library(dmrg-deps INTERFACE)
endif ()


# Setup dependencies
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)

# Install dependencies that are not in conan.
include(cmake/InstallPackage.cmake)
install_package(primme MODULE)
if (DMRG_ENABLE_TBLIS)
    install_package(tblis MODULE)
    target_link_libraries(dmrg-deps INTERFACE tblis::tblis)
endif ()


# Configure Eigen
if (TARGET Eigen3::Eigen)
    target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    if (TARGET BLAS::BLAS)
        target_link_libraries(Eigen3::Eigen INTERFACE BLAS::BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE_STRICT)
        if (TARGET mkl::mkl)
            message(STATUS "Eigen will use MKL")
            target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
        else ()
            message(STATUS "Eigen will use OpenBLAS")
        endif ()
    endif ()
else ()
    message(FATAL_ERROR "Target not defined: Eigen3::Eigen")
endif ()