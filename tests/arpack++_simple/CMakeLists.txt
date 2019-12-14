cmake_minimum_required(VERSION 3.10)
project(test-arpack++_simple)

add_executable(test-arpack++_simple arpack++_simple.cpp)
set(CMAKE_VERBOSE_MAKEFILE OFF)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-arpack++_simple PRIVATE -static)                                             ### Static linkage
endif()


target_link_libraries(test-arpack++_simple PRIVATE arpack++)
target_link_libraries(test-arpack++_simple PRIVATE project-settings)

add_dependencies(test-arpack++_simple arpack++)
add_test(NAME test-arpack++_simple COMMAND test-arpack++_simple)

