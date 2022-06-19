cmake_minimum_required(VERSION 3.15)
project(benchmarks)

file(GLOB BENCHMARKS "*.cpp")

set(dmrg-bench-tensoridentity-src)
set(dmrg-bench-swapgates-src
        #        $<TARGET_OBJECTS:dmrg-tid-o>
        #        $<TARGET_OBJECTS:dmrg-math-tenx>
        #        $<TARGET_OBJECTS:dmrg-math-svd>
        #        $<TARGET_OBJECTS:dmrg-tensors-mps>
        #        $<TARGET_OBJECTS:dmrg-tensors-state>
        #        $<TARGET_OBJECTS:dmrg-tools-common-contraction>
        #        $<TARGET_OBJECTS:dmrg-tools-finite-multisite>
        #        $<TARGET_OBJECTS:dmrg-tools-finite-measure>
        )
set(dmrg-bench-swapgates-lnk
        dmrg-tid
        dmrg-config
        dmrg-qm-gate
        dmrg-math-tenx
        dmrg-tensors-state-finite
        dmrg-tools-finite-mps
        dmrg-tools-finite-h5
        dmrg-tools-common-h5
        )


foreach (bench ${BENCHMARKS})
    get_filename_component(bench_src ${bench} NAME)
    get_filename_component(bench_nwe ${bench} NAME_WE)
    set(bench_tgt dmrg-bench-${bench_nwe})
    add_executable(${bench_tgt} ${bench_src} ${${bench_tgt}-src})
    target_link_libraries(${bench_tgt} PRIVATE ${${bench_tgt}-lnk})
    set_target_properties(${bench_tgt} PROPERTIES LINK_WHAT_YOU_USE TRUE)
    target_include_directories(${bench_tgt} PRIVATE ${CMAKE_SOURCE_DIR}/source)
    target_compile_definitions(${bench_tgt} PRIVATE BENCH_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
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

