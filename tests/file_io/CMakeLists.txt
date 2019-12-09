cmake_minimum_required(VERSION 3.10)
project(test-file_io)

add_executable(test-file_io file_io.cpp
        ../../source/general/nmspc_tensor_extra.h
        )
set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(test-file_io PUBLIC ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-file_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-file_io PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(test-file_io PRIVATE h5pp::h5pp h5pp::deps h5pp::flags)

target_link_libraries(test-file_io PRIVATE project-settings)
add_test(NAME test-file_io COMMAND test-file_io)
add_dependencies(test-file_io h5pp::h5pp)