add_source_target(xdmrg++-qm-time
        TARGET_SOURCES
        time.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-math-tenx
        xdmrg++-math-linalg
        xdmrg++-tools-common-log
        )