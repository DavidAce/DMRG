add_source_target(dmrg-tensors-model-finite
        TARGET_SOURCES
        ModelFinite.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-qm-spin
        dmrg-math-svd
        dmrg-tensors-mpo
        dmrg-tools-finite-multisite
        dmrg-tools-common-log
        )
add_source_target(dmrg-tensors-model-infinite
        TARGET_SOURCES
        ModelInfinite.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-config
        dmrg-math-svd
        dmrg-math-tenx
        dmrg-tensors-mpo

        )

add_library(dmrg-tensors-model INTERFACE)
target_link_libraries(dmrg-tensors-model INTERFACE dmrg-tensors-model-finite dmrg-tensors-model-infinite)