cmake_minimum_required(VERSION 3.18)
project(fill VERSION 1.0.0 DESCRIPTION "Add missing data to hdf5 files")

include(${CMAKE_SOURCE_DIR}/cmake/environment.cmake)

# Create an executable
add_executable(fill
               source/main.cpp
               source/parse.cpp
               source/mpi/mpi-tools.cpp
               )
target_include_directories(fill PRIVATE source)

target_link_libraries(fill PRIVATE
                      dmrg-debug
                      dmrg-stacktrace
                      dmrg-tools-common-log
                      dmrg-tools-finite-h5
                      dmrg-tools-finite-measure
                      dmrg-tensors-mpo-h5tb # To read the hamiltonian table
                      dmrg-tid)
target_link_libraries(fill PRIVATE dmrg-flags)
target_include_directories(fill PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_compile_definitions(fill PRIVATE CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
find_package(CLI11 REQUIRED)
find_package(h5pp REQUIRED)
find_package(Backward REQUIRED)
find_package(MPI COMPONENTS CXX REQUIRED)
target_link_libraries(fill PUBLIC CLI11::CLI11)
target_link_libraries(fill PUBLIC h5pp::h5pp)
target_link_libraries(fill PUBLIC Backward::Backward)
target_link_libraries(fill PRIVATE MPI::MPI_CXX)

# Print summary of CMake configuration
include(${CMAKE_SOURCE_DIR}/cmake/PrintTargetInfo.cmake)
print_project_summary(fill)



