cmake_minimum_required(VERSION 3.10)
project(test-file_io)

add_executable(${PROJECT_NAME} file_io.cpp
        ../../source/general/nmspc_tensor_extra.h
        )

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
