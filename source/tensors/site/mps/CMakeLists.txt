add_source_target(dmrg-tensors-mps
        TARGET_SOURCES
        MpsSite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-debug
        dmrg-math-linalg
        dmrg-math-tenx
        dmrg-math-hash
        dmrg-tools-common-contraction
        dmrg-tools-common-log
        )