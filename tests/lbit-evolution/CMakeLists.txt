cmake_minimum_required(VERSION 3.15)
project(dmrg-test-lbit-evolution)
enable_testing()
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
#                      dmrg-tid
#                      dmrg-math-svd
#                      dmrg-math-rnd
#                      dmrg-math-tenx
#                      dmrg-tools-common-log
#                      dmrg-tools-common-contraction
                      dmrg-tools-common-plot
                      dmrg-algorithms-flbit
                      )
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_precompiled_headers(${PROJECT_NAME})

if(DMRG_ENABLE_COVERAGE)
    target_compile_options(${PROJECT_NAME} PRIVATE --coverage)
    target_link_options(${PROJECT_NAME} PRIVATE --coverage)
endif()

add_test(NAME ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
               COMMAND ${PROJECT_NAME} -t 1 -c ${PROJECT_SOURCE_DIR}/input-flbit.cfg )

