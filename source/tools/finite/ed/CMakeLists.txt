add_source_target(dmrg-tools-finite-ed
        TARGET_SOURCES
        ed.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-algorithms-status
        dmrg-math-eig
        dmrg-math-tenx
        dmrg-tensors-finite
        dmrg-tools-finite-measure
        dmrg-tools-finite-opt
        dmrg-tools-common-split
        )

