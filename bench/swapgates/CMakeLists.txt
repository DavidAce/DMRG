cmake_minimum_required(VERSION 3.15)
project(dmrg-bench-swapgates)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-tid
                      dmrg-config
                      dmrg-debug
                      dmrg-qm-gate
                      dmrg-math-tenx
                      dmrg-tensors-state-finite
                      dmrg-tools-finite-mps
                      dmrg-tools-finite-h5
                      dmrg-tools-common-h5
                      dmrg-tools-common-prof)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_WHAT_YOU_USE TRUE)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_link_precompiled_headers(${PROJECT_NAME})
target_compile_definitions(${PROJECT_NAME} PRIVATE BENCH_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

