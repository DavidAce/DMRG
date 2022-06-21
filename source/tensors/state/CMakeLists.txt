add_source_target(dmrg-tensors-state-finite
        TARGET_SOURCES
        StateFinite.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-debug
        dmrg-tensors-mps
        dmrg-math-tenx
        dmrg-tools-finite-measure
        dmrg-tools-finite-multisite
        dmrg-tools-common-contraction
        dmrg-tools-common-log
        )
add_source_target(dmrg-tensors-state-infinite
        TARGET_SOURCES
        StateInfinite.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-debug
        dmrg-math-tenx
        dmrg-tensors-mps
        dmrg-tools-common-contraction
        dmrg-tools-infinite-measure
        dmrg-tools-infinite-mps
        dmrg-tools-common-log
        dmrg-tools-common-split
        dmrg-tools-common-views
        )

add_library(dmrg-tensors-state INTERFACE)
target_link_libraries(dmrg-tensors-state INTERFACE
        dmrg-tensors-state-finite
        dmrg-tensors-state-infinite)