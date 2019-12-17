cmake_minimum_required(VERSION 3.10)
project(test-ceres_pthread CXX)
add_executable(test-ceres_pthread ceres_pthread.cpp base.cpp rosenbrock.cpp)

set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(test-ceres_pthread PUBLIC ${CMAKE_SOURCE_DIR}/source)

if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-ceres_pthread PRIVATE -static)                                             ### Static linkage
endif()

if(CONAN_TARGETS)
    target_link_libraries(test-ceres_pthread PRIVATE CONAN_PKG::ceres-solver)
else()
    target_link_libraries(test-ceres_pthread PRIVATE ceres)
    target_link_libraries(test-ceres_pthread PRIVATE Eigen3::Eigen ) # Put it last in case Eigen wants to use blas
endif()
target_compile_features(test-ceres_pthread  PRIVATE cxx_std_17)

target_link_libraries(test-ceres_pthread PRIVATE project-settings)

add_test(NAME test-ceres_pthread COMMAND test-ceres_pthread)
