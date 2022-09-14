cmake_minimum_required(VERSION 3.15)
project(dmrg-bench-schmidt)
enable_testing()
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-tid
                      dmrg-math-svd
                      dmrg-math-tenx
                      dmrg-math-hash
                      dmrg-tensors-mps
                      dmrg-tools-common-contraction
                      dmrg-tools-common-log
                      dmrg-tools-common-split)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_link_precompiled_headers(${PROJECT_NAME})
target_compile_definitions(${PROJECT_NAME} PRIVATE BENCH_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

