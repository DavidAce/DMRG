cmake_minimum_required(VERSION 3.15)
project(xdmrg++-bench-expval)
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/bench)
target_link_libraries(${PROJECT_NAME} PRIVATE xdmrg++-deps xdmrg++-flags)
target_link_libraries(${PROJECT_NAME} PRIVATE
                      xdmrg++-config
                      xdmrg++-tools-common-log
                      xdmrg++-tools-common-contraction
                      )
if (DMRG_ENABLE_TBLIS)
    target_compile_definitions(${PROJECT_NAME} PRIVATE DMRG_ENABLE_TBLIS)
endif ()

target_link_precompiled_headers(${PROJECT_NAME})

