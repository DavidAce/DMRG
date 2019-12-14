cmake_minimum_required(VERSION 3.10)
project(test-eigen_openmp)
add_executable(test-eigen_openmp eigen_openmp.cpp)

set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(test-eigen_openmp PUBLIC ${CMAKE_SOURCE_DIR}/source)

if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-eigen_openmp PRIVATE -static)                                             ### Static linkage
endif()

if(CONAN_TARGETS)
    target_link_libraries(test-eigen_openmp PRIVATE CONAN_PKG::Eigen3)
else()
    target_link_libraries(test-eigen_openmp PRIVATE Eigen3::Eigen )
endif()

target_link_libraries(test-eigen_openmp PRIVATE blas)
target_link_libraries(test-eigen_openmp PRIVATE project-settings)

add_test(NAME test-eigen_openmp COMMAND test-eigen_openmp)
