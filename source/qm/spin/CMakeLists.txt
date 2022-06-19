add_source_target(dmrg-qm-spin
        TARGET_SOURCES
        spin.cpp
        half.cpp
        one.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        )