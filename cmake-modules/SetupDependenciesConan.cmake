
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


##################################################################
### Install conan-modules/conanfile.txt dependencies          ###
### This uses conan to get spdlog/eigen3/h5pp/ceres           ###
###    ceres-solver/2.0.0@davidace/development                ###
###    h5pp/0.4.7@davidace/stable                             ###
###    eigen/3.3.7@davidace/patched                           ###
##################################################################

find_program (
        CONAN_COMMAND
        conan
        HINTS $ENV{CONAN_PREFIX}
        PATHS $ENV{CONAN_PREFIX}/envs/dmrg $ENV{HOME}/anaconda3/envs/dmrg
        PATH_SUFFIXES bin
)
include(cmake-modules/conan/conan.cmake)
conan_add_remote(NAME conan-dmrg INDEX 1
        URL https://api.bintray.com/conan/davidace/conan-dmrg)
conan_cmake_run(CONANFILE cmake-modules/conan/conanfile.txt
        CONAN_COMMAND ${CONAN_COMMAND}
        SETTINGS compiler.cppstd=17
        SETTINGS compiler.libcxx=libstdc++11
        BASIC_SETUP CMAKE_TARGETS
        BUILD missing)



##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(project-settings INTERFACE CONAN_PKG::ceres-solver)
target_link_libraries(project-settings INTERFACE CONAN_PKG::h5pp)
target_link_libraries(project-settings INTERFACE CONAN_PKG::Eigen3)
target_link_libraries(project-settings INTERFACE arpack++)
if(TARGET OpenMP)
    target_link_libraries(project-settings INTERFACE OpenMP)
else()
    message(WARNING "We need OpenMP when using conan libraries!")
    target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
endif()


include(cmake-modules/PrintTargetInfo.cmake)
print_target_info(arpack++)
print_target_info(arpack)
print_target_info(lapack)
print_target_info(blas)
print_target_info(lapacke)
print_target_info(mkl)
print_target_info(OpenBLAS)
print_target_info(gfortran)
# Conan targets
foreach(tgt ${CONAN_TARGETS})
    if(NOT "${tgt}" MATCHES "Eigen|h5pp|ceres|arpack")
        print_target_info(${tgt})
    endif()
endforeach()
print_target_info(project-settings)
