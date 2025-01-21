add_source_target(xdmrg++-tools-infinite-measure
        TARGET_SOURCES
        measure.cpp
        measure_ham.cpp
        measure_mom.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-math-tenx
        xdmrg++-algorithms-status
        xdmrg++-tensors-infinite
        xdmrg++-tools-common-contraction
        xdmrg++-tools-common-log
        )

