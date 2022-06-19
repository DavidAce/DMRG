add_source_target(dmrg-tensors-edges-finite
        TARGET_SOURCES
        EdgesFinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-config
        dmrg-tensors-env
        dmrg-tools-common-log
        )
add_source_target(dmrg-tensors-edges-infinite
        TARGET_SOURCES
        EdgesInfinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tensors-env
        dmrg-tools-common-log
        )

add_library(dmrg-tensors-edges INTERFACE)
target_link_libraries(dmrg-tensors-edges INTERFACE dmrg-tensors-edges-finite dmrg-tensors-edges-infinite)