add_subdirectory(site)
add_subdirectory(model)
add_subdirectory(edges)
add_subdirectory(state)

add_source_target(dmrg-tensors-finite
        TARGET_SOURCES
        TensorsFinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-math-tenx
        dmrg-tensors-mps
        dmrg-tensors-mpo
        dmrg-tensors-env
        dmrg-tensors-state-finite
        dmrg-tensors-model-finite
        dmrg-tensors-edges-finite
        dmrg-tools-finite-env
        dmrg-tools-finite-measure
        dmrg-tools-finite-mps
        dmrg-tools-finite-multisite
        dmrg-tools-finite-ops
        dmrg-tools-common-log
        )
add_source_target(dmrg-tensors-infinite
        TARGET_SOURCES
        TensorsInfinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tensors-mps
        dmrg-tensors-mpo
        dmrg-tensors-env
        dmrg-tensors-state-infinite
        dmrg-tensors-model-infinite
        dmrg-tensors-edges-infinite
        dmrg-tools-infinite-env
        dmrg-tools-infinite-measure
        dmrg-tools-infinite-mps
        dmrg-tools-common-log
        )

add_library(dmrg-tensors INTERFACE)
target_link_libraries(dmrg-tensors INTERFACE dmrg-tensors-finite dmrg-tensors-infinite)
