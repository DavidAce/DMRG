add_source_target(dmrg-qm-h5tb
        TARGET_SOURCES
        h5tb.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-debug
        dmrg-tools-common-log
        )