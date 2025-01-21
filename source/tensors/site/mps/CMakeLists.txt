add_source_target(xdmrg++-tensors-mps
        TARGET_SOURCES
        MpsSite.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-debug
        xdmrg++-math-linalg
        xdmrg++-math-tenx
        xdmrg++-math-hash
        xdmrg++-tools-common-contraction
        xdmrg++-tools-common-log
        )