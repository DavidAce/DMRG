cmake_minimum_required(VERSION 3.15)
project(dmrg-bench-rsvd)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      dmrg-tid
                      dmrg-math-svd
                      dmrg-math-tenx
                      dmrg-tools-common-log)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_WHAT_YOU_USE TRUE)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_link_libraries(${PROJECT_NAME} PRIVATE dmrg-deps dmrg-flags)
target_link_precompiled_headers(${PROJECT_NAME})
target_compile_definitions(${PROJECT_NAME} PRIVATE BENCH_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

