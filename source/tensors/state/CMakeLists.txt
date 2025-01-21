add_source_target(xdmrg++-tensors-state-finite
        TARGET_SOURCES
        StateFinite.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-debug
        xdmrg++-tensors-mps
        xdmrg++-math-tenx
        xdmrg++-tools-finite-measure
        xdmrg++-tools-finite-multisite
        xdmrg++-tools-common-contraction
        xdmrg++-tools-common-log
        )
add_source_target(xdmrg++-tensors-state-infinite
        TARGET_SOURCES
        StateInfinite.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-debug
        xdmrg++-math-tenx
        xdmrg++-tensors-mps
        xdmrg++-tools-common-contraction
        xdmrg++-tools-infinite-measure
        xdmrg++-tools-infinite-mps
        xdmrg++-tools-common-log
        xdmrg++-tools-common-split
        xdmrg++-tools-common-views
        )

add_library(xdmrg++-tensors-state INTERFACE)
target_link_libraries(xdmrg++-tensors-state INTERFACE
        xdmrg++-tensors-state-finite
        xdmrg++-tensors-state-infinite)