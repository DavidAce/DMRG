add_source_target(xdmrg++-tools-finite-env
        TARGET_SOURCES
        env.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-config
        xdmrg++-math-svd
        xdmrg++-math-tenx
        xdmrg++-tensors-finite
        xdmrg++-tools-common-log
        xdmrg++-tools-common-contraction
        )

