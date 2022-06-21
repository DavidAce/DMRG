add_source_target(dmrg-tid
        TARGET_SOURCES
        tid.cpp
        ur.cpp
        token.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )