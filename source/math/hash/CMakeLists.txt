add_source_target(xdmrg++-math-hash
        TARGET_SOURCES
        hash.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        )
