add_library(dmrg-algorithms-o OBJECT)

target_sources(dmrg-algorithms-o PRIVATE
        AlgorithmBase.cpp
        AlgorithmFinite.cpp
        AlgorithmInfinite.cpp
        AlgorithmLauncher.cpp
        idmrg.cpp
        fdmrg.cpp
        xdmrg.cpp
        itebd.cpp
        flbit.cpp
        )

target_link_libraries(dmrg-algorithms-o PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-algorithms-o REUSE_FROM dmrg-pch)
endif ()

add_library(dmrg-algorithms INTERFACE)
target_link_libraries(dmrg-algorithms INTERFACE $<TARGET_OBJECTS:dmrg-algorithms-o>)
target_link_libraries(dmrg-algorithms INTERFACE dmrg-deps dmrg-flags)
