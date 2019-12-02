
cmake_minimum_required(VERSION 3.10)
#project(arpack++_simple)

add_executable(TEST-arpack++_simple arpack++_simple.cpp)
set(CMAKE_VERBOSE_MAKEFILE OFF)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-arpack++_simple PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(TEST-arpack++_simple PRIVATE arpack++)
target_link_libraries(TEST-arpack++_simple PRIVATE project-settings)

add_dependencies(TEST-arpack++_simple arpack++)
add_test(NAME arpack++_simple_test COMMAND TEST-arpack++_simple)

