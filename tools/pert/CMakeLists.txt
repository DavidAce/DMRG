cmake_minimum_required(VERSION 3.18)
project(pert VERSION 1.0.0 DESCRIPTION "Calculates the variance of non-interacting eigenstates with a small interacting perturbation")

include(${CMAKE_SOURCE_DIR}/cmake/environment.cmake)

# Create an executable
add_executable(pert
               source/main.cpp
               )
target_include_directories(pert PRIVATE source)

target_link_libraries(pert PRIVATE
                      dmrg-debug
                      dmrg-stacktrace
                      dmrg-tools-common-log
                      dmrg-tools-finite-h5
                      dmrg-tools-finite-measure
                      dmrg-tensors-mpo-h5tb # To read the hamiltonian table
                      dmrg-tid)
target_link_libraries(pert PRIVATE dmrg-flags)
target_include_directories(pert PRIVATE ${CMAKE_SOURCE_DIR}/source)

find_package(MPI COMPONENTS CXX REQUIRED)
find_package(CLI11 REQUIRED)
find_package(h5pp REQUIRED)
find_package(Backward REQUIRED)
target_link_libraries(pert PUBLIC CLI11::CLI11)
target_link_libraries(pert PUBLIC h5pp::h5pp)
target_link_libraries(pert PUBLIC Backward::Backward)

# Print summary of CMake configuration
include(${CMAKE_SOURCE_DIR}/cmake/PrintTargetInfo.cmake)
print_project_summary(pert)



