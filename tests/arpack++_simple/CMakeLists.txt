cmake_minimum_required(VERSION 3.10)
project(test-arpack++_simple)

add_executable(test-arpack++_simple arpack++_simple.cpp)
target_link_libraries(test-arpack++_simple PRIVATE project-settings)
add_test(NAME test-arpack++_simple COMMAND test-arpack++_simple)

