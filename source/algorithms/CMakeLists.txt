add_library(dmrg-algorithms-o OBJECT)

target_sources(dmrg-algorithms-o PRIVATE
        AlgorithmBase.cpp
        AlgorithmFinite.cpp
        AlgorithmInfinite.cpp
        AlgorithmLauncher.cpp
        AlgorithmStatus.cpp
        idmrg.cpp
        fdmrg.cpp
        xdmrg.cpp
        itebd.cpp
        flbit.cpp
        )

target_link_libraries(dmrg-algorithms-o PUBLIC deps flags)
target_link_precompiled_headers(dmrg-algorithms-o)
add_library(dmrg-algorithms INTERFACE)
target_link_libraries(dmrg-algorithms INTERFACE $<TARGET_OBJECTS:dmrg-algorithms-o>)
target_link_libraries(dmrg-algorithms INTERFACE deps flags)
