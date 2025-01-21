add_source_target(xdmrg++-tools-finite-measure
        TARGET_SOURCES
        measure.cpp
        measure-parity-components.cpp
        entropy.cpp
        information.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-qm-spin
        xdmrg++-qm-mpo
        xdmrg++-tensors-finite
        xdmrg++-math-eig
        xdmrg++-math-tenx
        xdmrg++-math-linalg
        xdmrg++-tools-finite-mps
        xdmrg++-tools-common-log
        xdmrg++-tools-common-contraction

        )

