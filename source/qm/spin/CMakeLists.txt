add_source_target(dmrg-qm-spin
        TARGET_SOURCES
        spin.cpp
        half.cpp
        one.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )