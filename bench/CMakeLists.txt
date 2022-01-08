cmake_minimum_required(VERSION 3.15)
project(benchmarks)
file(GLOB BENCHMARKS "*.cpp")
foreach (bench ${BENCHMARKS})
    get_filename_component(bench_src ${bench} NAME)
    get_filename_component(bench_nwe ${bench} NAME_WE)
    set(bench_tgt dmrg-${bench_nwe})
    add_executable(${bench_tgt} ${bench_src} ${bench_inc})
    target_link_libraries(${bench_tgt} PUBLIC ${bench_lnk})
    set_target_properties(${bench_tgt} PROPERTIES LINK_WHAT_YOU_USE TRUE)
    target_include_directories(${bench_tgt} PRIVATE ${CMAKE_SOURCE_DIR}/source)
    target_compile_definitions(${bench_tgt} PRIVATE BENCH_MATRIX_DIR="${PROJECT_SOURCE_DIR}/benchmatrices")
    add_test(NAME ${bench_tgt} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${bench_tgt})
    if (MSVC)
        set_target_properties(${bench_tgt} PROPERTIES LINK_FLAGS "/ignore:4099")
    endif ()
    list(APPEND bench_tgt_list ${bench_tgt})
    target_link_libraries(${bench_tgt} PUBLIC deps flags)
    # Speed up compilation with precompiled headers
    target_link_precompiled_headers(${bench_tgt})
endforeach ()


add_custom_target(dmrg-bench-all
        COMMENT "Building benchmarks"
        DEPENDS ${bench_tgt_list}
        )
set(dmrg-bench-list ${bench_tgt_list} PARENT_SCOPE)

