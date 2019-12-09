
cmake_minimum_required(VERSION 3.10)
project(test-arpack++_mps)

enable_language(Fortran)
set(CMAKE_VERBOSE_MAKEFILE OFF)
add_executable(test-arpack++_mps arpack++_mps.cpp
                            ../../source/math/class_eigsolver.cpp
                            ../../source/math/arpack_extra/arpackpp_solver.cpp
                            ../../source/math/arpack_extra/matrix_product_stl.cpp
                            ../../source/general/class_tic_toc.cpp)
target_include_directories(test-arpack++_mps PUBLIC ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-arpack++_mps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-arpack++_mps PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(test-arpack++_mps PRIVATE project-settings)
target_link_libraries(test-arpack++_mps PRIVATE arpack++ Eigen3::Eigen spdlog::spdlog)

add_test(NAME test-arpack++_mps COMMAND test-arpack++_mps)
add_dependencies(test-arpack++_mps arpack++ Eigen3::Eigen spdlog::spdlog)