cmake_minimum_required(VERSION 3.15)
project(dmrg-bench-tensoridentity)
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_link_precompiled_headers(${PROJECT_NAME})

