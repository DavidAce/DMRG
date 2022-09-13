cmake_minimum_required(VERSION 3.15)
project(dmrg-test-map-shuffle)
enable_testing()
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-tid
                      dmrg-math-svd
                      dmrg-math-rnd
                      dmrg-math-tenx
                      dmrg-tools-common-log
                      dmrg-tools-common-contraction
                      )
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_WHAT_YOU_USE TRUE)
target_link_precompiled_headers(${PROJECT_NAME})

if(DMRG_ENABLE_COVERAGE)
    target_compile_options(${PROJECT_NAME} PRIVATE --coverage)
    target_link_options(${PROJECT_NAME} PRIVATE --coverage)
endif()

add_test(NAME ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${PROJECT_NAME})

