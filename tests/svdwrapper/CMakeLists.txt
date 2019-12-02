
cmake_minimum_required(VERSION 3.10)
#project(svdwrapper)
add_executable(TEST-svdwrapper svdwrapper.cpp
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper.h
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper.cpp
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper_lapacke.cpp)

set(CMAKE_VERBOSE_MAKEFILE OFF)

target_include_directories(TEST-svdwrapper PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-svdwrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(TEST-svdwrapper PRIVATE TEST_DIR="${PROJECT_SOURCE_DIR}/tests/svdwrapper")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-svdwrapper PRIVATE -static)                                             ### Static linkage
endif()
target_link_libraries(TEST-svdwrapper PRIVATE h5pp::h5pp h5pp::deps h5pp::flags lapacke)
target_link_libraries(TEST-svdwrapper PRIVATE project-settings)
add_test(NAME svdwrapper_test COMMAND TEST-svdwrapper)

add_dependencies(TEST-svdwrapper h5pp::h5pp h5pp::deps h5pp::flags gfortran )