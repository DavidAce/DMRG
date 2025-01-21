add_source_target(xdmrg++-tools-infinite-env
        TARGET_SOURCES
        env.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tensors-infinite
        )

