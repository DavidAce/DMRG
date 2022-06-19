add_source_target(dmrg-tools-infinite-env
        TARGET_SOURCES
        env.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tensors-infinite
        )

