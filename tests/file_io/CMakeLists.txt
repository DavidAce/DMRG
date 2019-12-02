cmake_minimum_required(VERSION 3.10)
#project(file_io)

add_executable(TEST-file_io file_io.cpp
        ../../source/general/nmspc_tensor_extra.h
        )
set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(TEST-file_io PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-file_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-file_io PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(TEST-file_io PRIVATE h5pp::h5pp h5pp::deps h5pp::flags)

target_link_libraries(TEST-file_io PRIVATE project-settings)
add_test(NAME file_io_test COMMAND TEST-file_io)
add_dependencies(TEST-file_io h5pp::h5pp)