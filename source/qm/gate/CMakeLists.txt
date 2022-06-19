add_source_target(dmrg-qm-gate
        TARGET_SOURCES
        gate.cpp
        swap.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-math-tenx
        dmrg-tools-common-log
        )