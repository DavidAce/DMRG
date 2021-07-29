add_library(dmrg-tools-finite-ops OBJECT ops.cpp)
add_library(dmrg-tools-finite-print OBJECT print.cpp)

target_link_libraries(dmrg-tools-finite-ops PUBLIC dmrg-deps dmrg-flags)
target_link_libraries(dmrg-tools-finite-print PUBLIC dmrg-deps dmrg-flags)

if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-tools-finite-ops REUSE_FROM dmrg-pch)
    target_precompile_headers(dmrg-tools-finite-print REUSE_FROM dmrg-pch)
endif ()

add_subdirectory(h5)
add_subdirectory(measure)
add_subdirectory(multisite)
add_subdirectory(env)
add_subdirectory(mps)
add_subdirectory(opt)

add_library(dmrg-tools-finite INTERFACE)
target_link_libraries(dmrg-tools-finite INTERFACE
        $<TARGET_OBJECTS:dmrg-tools-finite-env>
        $<TARGET_OBJECTS:dmrg-tools-finite-mps>
        $<TARGET_OBJECTS:dmrg-tools-finite-multisite>
        $<TARGET_OBJECTS:dmrg-tools-finite-ops>
        $<TARGET_OBJECTS:dmrg-tools-finite-opt>
        $<TARGET_OBJECTS:dmrg-tools-finite-print>
        $<TARGET_OBJECTS:dmrg-tools-finite-h5>
        $<TARGET_OBJECTS:dmrg-tools-finite-measure>
        )
target_link_libraries(dmrg-tools-finite INTERFACE dmrg-deps dmrg-flags)
