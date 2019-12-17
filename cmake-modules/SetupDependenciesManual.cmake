
##############################################################################
###  Optional OpenMP support                                               ###
###  Note that Clang has some  trouble with static openmp and that         ###
###  and that static openmp is not recommended. This tries to enable       ###
###  static openmp anyway because I find it useful. Installing             ###
###  libiomp5 might help for shared linking.                               ###
##############################################################################
include(cmake-modules/FindOpenMPLibrary.cmake)
find_package_openmp()



# These packages are not in conan yet
include(cmake-modules/FindGFortran.cmake)
if(ENABLE_MKL)
    include(cmake-modules/Find_dont_install_INTELMKL.cmake)    # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
endif()
if(NOT TARGET blas)
    include(cmake-modules/Fetch_OpenBLAS.cmake)                 # If MKL is not on openblas will be used instead. Includes lapack.
endif()

include(cmake-modules/FindLapacke.cmake)                        # Lapacke needed by arpack++
include(cmake-modules/Fetch_arpack-ng.cmake)                    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake-modules/Fetch_arpack++.cmake)                     # LC++ frontend for arpack-ng
include(cmake-modules/Fetch_Eigen3.cmake)                       # Eigen3 numerical library (needed by ceres and h5pp)
include(cmake-modules/Fetch_h5pp.cmake)                         # GSL for computing some numerical integrals.
include(cmake-modules/Fetch_gflags.cmake)                       # Google Flags library needed by ceres-solver
include(cmake-modules/Fetch_glog.cmake)                         # Google logging library needed by ceres-solver
include(cmake-modules/Fetch_ceres-solver.cmake)                 # ceres-solver (for L-BFGS routine)


##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(project-settings INTERFACE ceres)
target_link_libraries(project-settings INTERFACE h5pp::h5pp h5pp::deps h5pp::flags)
target_link_libraries(project-settings INTERFACE arpack++)
target_link_libraries(project-settings INTERFACE Eigen3::Eigen) # Put it last in case Eigen wants to use blas
if(TARGET OpenMP)
    target_link_libraries(project-settings INTERFACE OpenMP)
else()
    target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
endif()
target_link_libraries(project-settings INTERFACE -lrt -Wl,--whole-archive  -lpthread  -Wl,--no-whole-archive)



include(cmake-modules/PrintTargetInfo.cmake)
print_target_info(ceres)
print_target_info(gflags)
print_target_info(glog::glog)
print_target_info(h5pp::h5pp)
print_target_info(h5pp::deps)
print_target_info(h5pp::flags)
print_target_info(spdlog::spdlog)
print_target_info(hdf5::hdf5)
print_target_info(Eigen3::Eigen)
print_target_info(arpack++)
print_target_info(arpack)
print_target_info(lapack)
print_target_info(blas)
print_target_info(lapacke)
print_target_info(mkl)
print_target_info(OpenBLAS)
print_target_info(gfortran)
print_target_info(Threads::Threads)
print_target_info(OpenMP)
print_target_info(project-settings)



