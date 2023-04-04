
##############################################################################
###  Optional Intel MKL support. Uses OpenBLAS as fall-back                ###
##############################################################################
if(NOT DEFINED BLA_VENDOR AND NOT DEFINED $ENV{BLA_VENDOR})
    message(FATAL_ERROR  "BLA_VENDOR is undefined: https://cmake.org/cmake/help/latest/module/FindBLAS.html" )
endif()

if(BLA_VENDOR MATCHES Intel OR ENV{BLA_VENDOR} MATCHES Intel)
    find_package(MKL REQUIRED BYPASS_PROVIDER)
endif()

# Link all dependencies to dmrg-deps
if(NOT TARGET dmrg-deps)
    add_library(dmrg-deps INTERFACE)
endif()

find_package(OpenMP COMPONENTS C CXX REQUIRED BYPASS_PROVIDER)


# Setup dependencies
include(cmake/SetupDependenciesFind.cmake)
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)

# Install dependencies that are not in conan.
include(cmake/InstallPackage.cmake)

install_package(primme MODULE)
target_link_libraries(dmrg-deps INTERFACE primme::primme)

if(DMRG_ENABLE_TBLIS)
    install_package(tblis MODULE)
    target_link_libraries(dmrg-deps INTERFACE tblis::tblis)
endif()

# Configure Eigen
if(TARGET Eigen3::Eigen)
    target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_THREADS)
    if(TARGET MKL::MKL)
        message(STATUS "Eigen will use MKL")
        target_link_libraries(Eigen3::Eigen INTERFACE MKL::MKL)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
    elseif(TARGET BLAS::BLAS)
        message(STATUS "Eigen will use BLAS")
        target_link_libraries(Eigen3::Eigen INTERFACE BLAS::BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE)
    endif()
else()
    message(FATAL_ERROR "Target not defined: Eigen3::Eigen")
endif()

#set_target_properties(OpenMP::OpenMP_CXX PROPERTIES INTERFACE_LINK_LIBRARIES "${OpenMP_CXX_FLAGS}") # Use flag only
