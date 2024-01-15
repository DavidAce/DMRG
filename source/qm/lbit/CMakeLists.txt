add_source_target(dmrg-qm-lbit-h5tb
        TARGET_SOURCES
        h5tb.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )
add_source_target(dmrg-qm-lbit
        TARGET_SOURCES
        lbit.cpp
        cls.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-qm-lbit-h5tb
        dmrg-deps dmrg-flags
        )

