add_source_target(dmrg-tools-finite-env
        TARGET_SOURCES
        env.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-math-svd
        dmrg-math-tenx
        dmrg-tensors-finite
        dmrg-tools-common-log
        dmrg-tools-common-contraction
        )

