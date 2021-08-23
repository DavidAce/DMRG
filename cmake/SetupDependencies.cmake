
#################################################################
### Preempt Threads::Threads                                   ###
### It's looked for in dependencies, so we make it right       ###
### before it's done wrong, i.e. with pthread instead of       ###
### -lpthread.                                                 ###
### Here we specify the linking twice                          ###
### 1) As string to make sure -lpthread gets sandwiched by     ###
###    -Wl,--whole-archive.... -Wl,--no-whole-archive          ###
###    -Wl,--whole-archive.... -Wl,--no-whole-archive          ###
### 2) As usual to make sure that if somebody links            ###
###    Threads::Threads, then any repeated pthread appended    ###
###    to the end (the wrong order causes linking errors)      ###
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
    target_link_libraries(dmrg-flags INTERFACE OpenMP::OpenMP_CXX)
else()
    target_compile_options(dmrg-flags INTERFACE -Wno-unknown-pragmas)
endif()

add_library(dmrg-deps INTERFACE)
target_link_libraries(dmrg-deps INTERFACE
        cxxopts::cxxopts
        h5pp::h5pp
        arpack::arpack++
        primme::primme
        Ceres::ceres
        )
if(TARGET unwind::unwind)
    target_link_libraries(dmrg-deps INTERFACE unwind::unwind)
    target_compile_definitions(dmrg-deps INTERFACE DMRG_HAS_UNWIND=1)
endif()