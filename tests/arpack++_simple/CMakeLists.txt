cmake_minimum_required(VERSION 3.10)
project(test-arpack++_simple)

add_executable(${PROJECT_NAME} arpack++_simple.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

