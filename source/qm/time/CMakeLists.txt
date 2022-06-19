add_source_target(dmrg-qm-time
        TARGET_SOURCES
        time.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-math-tenx
        dmrg-math-linalg
        dmrg-tools-common-log
        )