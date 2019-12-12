
cmake_minimum_required(VERSION 3.10)
project(test-svdwrapper)
add_executable(test-svdwrapper svdwrapper.cpp
        ${CMAKE_SOURCE_DIR}/source/math/class_svd_wrapper.h
        ${CMAKE_SOURCE_DIR}/source/math/class_svd_wrapper.cpp
        ${CMAKE_SOURCE_DIR}/source/math/class_svd_wrapper_lapacke.cpp)

set(CMAKE_VERBOSE_MAKEFILE OFF)

target_include_directories(test-svdwrapper PUBLIC ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-svdwrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(test-svdwrapper PRIVATE TEST_DIR="${PROJECT_SOURCE_DIR}")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-svdwrapper PRIVATE -static)                                             ### Static linkage
endif()
target_link_libraries(test-svdwrapper PRIVATE Eigen3::Eigen h5pp::h5pp h5pp::deps h5pp::flags)
target_link_libraries(test-svdwrapper PRIVATE project-settings)
add_test(NAME test-svdwrapper COMMAND test-svdwrapper)
