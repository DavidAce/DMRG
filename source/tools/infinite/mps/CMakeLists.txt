add_source_target(xdmrg++-tools-infinite-mps
        TARGET_SOURCES
        mps.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-config
        xdmrg++-tensors-state-infinite
        xdmrg++-tools-common-split
        )

