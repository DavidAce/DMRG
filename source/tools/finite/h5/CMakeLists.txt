add_source_target(dmrg-tools-finite-h5
        TARGET_SOURCES
        load.cpp
        save.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-algorithms-status
        dmrg-math-tenx
        dmrg-tensors-finite
        dmrg-tools-finite-measure
        dmrg-tools-finite-ops
        dmrg-tools-finite-env
        dmrg-tools-finite-print
        dmrg-tools-common-h5
        dmrg-tools-common-log
        )

