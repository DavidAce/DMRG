cmake_minimum_required(VERSION 3.10)
project(test-eigen_openmp)
add_executable(test-eigen_openmp eigen_openmp.cpp)

set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(test-eigen_openmp PUBLIC ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(test-eigen_openmp PRIVATE project-settings)
add_test(NAME test-eigen_openmp COMMAND test-eigen_openmp)
