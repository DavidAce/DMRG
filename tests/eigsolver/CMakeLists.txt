
cmake_minimum_required(VERSION 3.10)
#project(eigsolver)
add_executable(TEST-eigsolver eigsolver.cpp
        ${PROJECT_SOURCE_DIR}/source/math/class_eigsolver.cpp
        ${PROJECT_SOURCE_DIR}/source/general/class_tic_toc.cpp)

set(CMAKE_VERBOSE_MAKEFILE OFF)

target_include_directories(TEST-eigsolver PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-eigsolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(TEST-eigsolver PRIVATE TEST_MATRIX_DIR="${PROJECT_SOURCE_DIR}/tests/eigsolver")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-eigsolver PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(TEST-eigsolver PRIVATE project-settings)
target_link_libraries(TEST-eigsolver PRIVATE h5pp::h5pp h5pp::deps h5pp::flags lapacke )
add_test(NAME eigsolver_test COMMAND TEST-eigsolver)

add_dependencies(TEST-eigsolver h5pp::h5pp h5pp::deps h5pp::flags gfortran )