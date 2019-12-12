
cmake_minimum_required(VERSION 3.10)
project(test-eigsolver)
add_executable(test-eigsolver eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/math/class_eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/general/class_tic_toc.cpp)

set(CMAKE_VERBOSE_MAKEFILE OFF)

target_include_directories(test-eigsolver PUBLIC ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-eigsolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(test-eigsolver PRIVATE TEST_MATRIX_DIR="${PROJECT_SOURCE_DIR}")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-eigsolver PRIVATE -static)                                             ### Static linkage
endif()
target_link_libraries(test-eigsolver PRIVATE blas lapack lapacke Eigen3::Eigen h5pp::h5pp h5pp::deps h5pp::flags)
target_link_libraries(test-eigsolver PRIVATE project-settings)
add_test(NAME test-eigsolver COMMAND test-eigsolver)
