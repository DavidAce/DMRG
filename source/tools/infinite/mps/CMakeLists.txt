add_source_target(dmrg-tools-infinite-mps
        TARGET_SOURCES
        mps.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-tensors-state-infinite
        dmrg-tools-common-split
        )

