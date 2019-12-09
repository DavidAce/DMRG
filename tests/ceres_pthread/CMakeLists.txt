cmake_minimum_required(VERSION 3.10)
project(test-ceres_pthread)
add_executable(test-ceres_pthread ceres_pthread.cpp)

set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(test-ceres_pthread PUBLIC ${CMAKE_SOURCE_DIR}/source)

if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (test-ceres_pthread PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(test-ceres_pthread PRIVATE ceres)
target_link_libraries(test-ceres_pthread PRIVATE project-settings)

add_test(NAME test-ceres_pthread COMMAND test-ceres_pthread)
