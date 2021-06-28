#include(cmake/SetupDependenciesFind.cmake)
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)


include(cmake/InstallPackage.cmake)
install_package(primme MODULE)

target_link_libraries(dmrg-eig PUBLIC primme::primme)

##################################################################
### Link all the things!                                       ###
##################################################################
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(dmrg-flags INTERFACE OpenMP::OpenMP_CXX)
else()
    target_compile_options(dmrg-flags INTERFACE -Wno-unknown-pragmas)
endif()
target_link_libraries(dmrg-deps INTERFACE h5pp::h5pp arpack::arpack++ Ceres::ceres primme::primme)
target_link_libraries(dmrg-main PUBLIC h5pp::h5pp)
target_link_libraries(dmrg-opt PUBLIC spdlog::spdlog Ceres::ceres)
target_link_libraries(dmrg-eig PUBLIC spdlog::spdlog Eigen3::Eigen)
target_link_libraries(dmrg-arp PUBLIC spdlog::spdlog arpack::arpack++ primme::primme)

if(TARGET unwind::unwind)
    target_link_libraries(dmrg-dbg PUBLIC unwind::unwind)
    target_compile_definitions(dmrg-dbg PUBLIC DMRG_HAS_UNWIND=1)
endif()