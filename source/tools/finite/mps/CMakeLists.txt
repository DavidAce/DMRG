add_source_target(dmrg-tools-finite-mps
        TARGET_SOURCES
        mps.cpp
        mps_init_entangled_state.cpp
        mps_init_product_state.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-math-svd
        dmrg-qm-gate
        dmrg-qm-mpo
        dmrg-qm-spin
        dmrg-tensors-state-finite
        dmrg-tensors-mps
        dmrg-tools-common-contraction
        dmrg-tools-common-log
        dmrg-tools-common-split
        dmrg-tools-finite-measure
        dmrg-tools-finite-ops
        )

