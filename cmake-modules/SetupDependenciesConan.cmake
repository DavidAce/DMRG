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
if(NOT TARGET blas::blas)
    include(cmake-modules/Fetch_OpenBLAS.cmake)                 # If MKL is not on openblas will be used instead. Includes lapack.
endif()
include(cmake-modules/FindLapacke.cmake)                        # Lapacke needed by arpack++
include(cmake-modules/Fetch_arpack-ng.cmake)                    # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake-modules/Fetch_arpack++.cmake)                     # LC++ frontend for arpack-ng


##################################################################
### Install conan-modules/conanfile.txt dependencies          ###
### This uses conan to get spdlog/eigen3/h5pp/ceres           ###
###    ceres-solver/2.0.0@davidace/development                ###
###    h5pp/1.1.0@davidace/stable                             ###
###    eigen/3.3.7@davidace/patched                           ###
##################################################################

find_program (
        CONAN_COMMAND
        conan
        HINTS ${CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONAN_PREFIX} $ENV{CONDA_PREFIX}
        PATHS $ENV{HOME}/anaconda3 $ENV{HOME}/miniconda $ENV{HOME}/.conda
        PATH_SUFFIXES bin envs/dmrg/bin
)


# Download automatically, you can also just copy the conan.cmake file
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.14/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake")
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)

conan_add_remote(NAME conan-dmrg INDEX 1
        URL https://api.bintray.com/conan/davidace/conan-dmrg)


if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    # Let it autodetect libcxx
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # There is no libcxx
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    list(APPEND conan_libcxx compiler.libcxx=libstdc++11)
endif()


conan_cmake_run(CONANFILE conanfile.txt
        CONAN_COMMAND ${CONAN_COMMAND}
        SETTINGS compiler.cppstd=17
        SETTINGS "${conan_libcxx}"
        BUILD_TYPE ${CMAKE_BUILD_TYPE}
        OUTPUT_QUIET
        BASIC_SETUP CMAKE_TARGETS
        BUILD missing)



##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(project-settings INTERFACE CONAN_PKG::ceres-solver)
target_link_libraries(project-settings INTERFACE CONAN_PKG::h5pp)
target_link_libraries(project-settings INTERFACE CONAN_PKG::Eigen3)
target_link_libraries(project-settings INTERFACE arpack::arpack++) # Last to use blas


if(TARGET openmp::openmp)
    target_link_libraries(project-settings INTERFACE openmp::openmp)
else()
    message(WARNING "We need OpenMP when using conan libraries!")
    target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
endif()

target_link_libraries(project-settings INTERFACE -Wl,--whole-archive pthread -Wl,--no-whole-archive -lrt -ldl )

include(cmake-modules/PrintTargetInfo.cmake)
print_target_info(arpack::arpack++)
print_target_info(arpack::arpack)
print_target_info(lapack::lapack)
print_target_info(blas::blas)
print_target_info(lapacke::lapacke)

# Conan targets
foreach(tgt ${CONAN_TARGETS})
    print_target_info(${tgt})
endforeach()
print_target_info(openmp::openmp)
print_target_info(mkl::mkl)
print_target_info(openblas::openblas)
print_target_info(gfortran::gfortran)
print_target_info(Threads::Threads)
print_target_info(project-settings)
