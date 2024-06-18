add_source_target(dmrg-tools-finite-measure
        TARGET_SOURCES
        measure.cpp
        measure-parity-components.cpp
        entropy.cpp
        information.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-qm-spin
        dmrg-qm-mpo
        dmrg-tensors-finite
        dmrg-math-eig
        dmrg-math-tenx
        dmrg-math-linalg
        dmrg-tools-finite-mps
        dmrg-tools-common-log
        dmrg-tools-common-contraction

        )

