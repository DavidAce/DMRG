
cmake_minimum_required(VERSION 3.10)
project(test-eigsolver)
add_executable(${PROJECT_NAME} eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/math/class_eigsolver.cpp
        ${CMAKE_SOURCE_DIR}/source/general/class_tic_toc.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(${PROJECT_NAME} PRIVATE TEST_MATRIX_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
