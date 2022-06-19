add_source_target(dmrg-qm-mpo
        TARGET_SOURCES
        mpo.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-math-rnd
        dmrg-math-tenx
        )
