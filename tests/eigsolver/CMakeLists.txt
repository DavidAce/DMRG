
cmake_minimum_required(VERSION 3.10)
project(test-eigsolver)
add_executable(test-eigsolver eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/math/class_eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/general/class_tic_toc.cpp)

target_include_directories(test-eigsolver PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-eigsolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(test-eigsolver PRIVATE TEST_MATRIX_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(test-eigsolver PRIVATE project-settings)
add_test(NAME test-eigsolver COMMAND test-eigsolver)
