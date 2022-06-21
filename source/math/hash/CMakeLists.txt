add_source_target(dmrg-math-hash
        TARGET_SOURCES
        hash.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        )
