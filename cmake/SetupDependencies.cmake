
##################################################################
### Preempt Threads::Threads                                   ###
### It's looked for in dependencies, so we make it right       ###
### before it's done wrong, i.e. with pthread instead of       ###
### -lpthread which causes link errors downstream with         ###
###    -Wl,--whole-archive.... -Wl,--no-whole-archive          ###
##################################################################
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Threads::Threads INTERFACE rt dl)


include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)

include(cmake/InstallPackage.cmake)
install_package(primme MODULE)


##################################################################
### Link all the things!                                       ###
##################################################################
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(flags INTERFACE OpenMP::OpenMP_CXX)
else()
    target_compile_options(flags INTERFACE -Wno-unknown-pragmas)
endif()

if(NOT TARGET deps)
    add_library(deps INTERFACE)
endif()
target_link_libraries(deps INTERFACE
        CLI11::CLI11
        h5pp::h5pp
        arpack::arpack++
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