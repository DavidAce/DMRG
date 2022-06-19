add_source_target(dmrg-tools-infinite-h5
        TARGET_SOURCES
        save.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-algorithms-status
        dmrg-tensors-infinite
        dmrg-tools-infinite-measure
        )

