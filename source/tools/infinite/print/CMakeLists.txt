add_source_target(xdmrg++-tools-infinite-print
        TARGET_SOURCES
        print.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tensors-infinite
        )

