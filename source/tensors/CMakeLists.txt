add_library(dmrg-tensors-tensors OBJECT)
target_sources(dmrg-tensors-tensors PRIVATE TensorsFinite.cpp TensorsInfinite.cpp)
target_link_libraries(dmrg-tensors-tensors PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-tensors)


add_subdirectory(site)
add_subdirectory(model)
add_subdirectory(edges)
add_subdirectory(state)

add_library(dmrg-tensors INTERFACE)
target_sources(dmrg-tensors INTERFACE
        $<TARGET_OBJECTS:dmrg-tensors-mps>
        $<TARGET_OBJECTS:dmrg-tensors-mpo>
        $<TARGET_OBJECTS:dmrg-tensors-env>
        $<TARGET_OBJECTS:dmrg-tensors-state>
        $<TARGET_OBJECTS:dmrg-tensors-model>
        $<TARGET_OBJECTS:dmrg-tensors-edges>
        $<TARGET_OBJECTS:dmrg-tensors-tensors>
        )
