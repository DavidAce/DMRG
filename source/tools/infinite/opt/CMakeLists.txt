add_source_target(dmrg-tools-infinite-opt
        TARGET_SOURCES
        opt.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-eig
        dmrg-math-tenx
        dmrg-tensors-infinite
        dmrg-tools-common-split
        )

