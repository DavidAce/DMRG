add_source_target(xdmrg++-qm-spin
        TARGET_SOURCES
        spin.cpp
        half.cpp
        one.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        )