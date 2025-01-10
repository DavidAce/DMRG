cmake_minimum_required(VERSION 3.15)
project(dmrg-bench-eig)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-tid
                      dmrg-math
                      dmrg-debug
                      dmrg-tools-common-log
                      dmrg-tools-common-contraction
                      )
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_compile_definitions(${PROJECT_NAME} PRIVATE BENCH_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
target_link_precompiled_headers(${PROJECT_NAME})

