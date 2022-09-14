cmake_minimum_required(VERSION 3.15)
project(dmrg-test-partial-trace)
enable_testing()
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-qm-spin
                      dmrg-math-rnd
                      dmrg-math-linalg
                      )
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_precompiled_headers(${PROJECT_NAME})

if(DMRG_ENABLE_COVERAGE)
    target_compile_options(${PROJECT_NAME} PRIVATE --coverage)
    target_link_options(${PROJECT_NAME} PRIVATE --coverage)
endif()

add_test(NAME ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${PROJECT_NAME})
