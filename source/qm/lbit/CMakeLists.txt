add_source_target(dmrg-qm-lbit
        TARGET_SOURCES
        lbit.cpp
        cls.cpp
        h5tb.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )