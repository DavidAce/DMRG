add_source_target(dmrg-tools-finite-multisite
        TARGET_SOURCES
        multisite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tensors-state-finite
        dmrg-tensors-model-finite
        dmrg-tools-finite-measure
        dmrg-tools-finite-mps
        dmrg-tools-common-log
        dmrg-tools-common-contraction
        )

